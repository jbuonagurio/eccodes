/*
 * Copyright 2005-2018 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
 * virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
 */

#include "grib_api_internal.h"

#if GRIB_PTHREADS
 static pthread_once_t once  = PTHREAD_ONCE_INIT;
 static pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
 static pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
 static void init() {
  pthread_mutexattr_t attr;
  pthread_mutexattr_init(&attr);
  pthread_mutexattr_settype(&attr,PTHREAD_MUTEX_RECURSIVE);
  pthread_mutex_init(&mutex1,&attr);
  pthread_mutex_init(&mutex2,&attr);
  pthread_mutexattr_destroy(&attr);
 }
#elif GRIB_OMP_THREADS
 static int once = 0;
 static omp_nest_lock_t mutex1;
 static omp_nest_lock_t mutex2;
 static void init()
 {
    GRIB_OMP_CRITICAL(lock_grib_io_c)
    {
        if (once == 0)
        {
            omp_init_nest_lock(&mutex1);
            omp_init_nest_lock(&mutex2);
            once = 1;
        }
    }
 }
#endif

#ifdef HAVE_CURL
 #include <curl/curl.h>
#endif

#define  GRIB 0x47524942
#define  BUDG 0x42554447
#define  DIAG 0x44494147
#define  TIDE 0x54494445
#define  BUFR 0x42554652
#define  HDF5 0x89484446
#define  WRAP 0x57524150

#define GRIB_API_READS_BUFR 1
#define GRIB_API_READS_HDF5 1
#define GRIB_API_READS_WRAP 1


typedef struct alloc_buffer {
    size_t size;
    void* buffer;
} alloc_buffer;

typedef size_t   (*readproc)(void*,void*,size_t,int*);
typedef int   (*seekproc)(void*,off_t);
typedef off_t   (*tellproc)(void*);
typedef void* (*allocproc)(void*,size_t*,int*);


typedef struct reader {
    void *read_data;
    readproc read;

    void *alloc_data;
    allocproc alloc;
    int headers_only;

    seekproc seek;
    seekproc seek_from_start;
    tellproc tell;
    off_t offset;

    size_t message_size;

} reader;

static int read_the_rest(reader* r,size_t message_length,unsigned char* tmp, int already_read, int check7777)
{
    int err = GRIB_SUCCESS;
    size_t buffer_size;
    size_t rest;
    unsigned char* buffer;

    if (message_length==0)
        return GRIB_BUFFER_TOO_SMALL;

    buffer_size = message_length;
    rest=message_length-already_read;
    r->message_size=message_length;
    buffer = (unsigned char*)r->alloc(r->alloc_data,&buffer_size,&err);
    if(err) return err;

    if (buffer == NULL || (buffer_size < message_length)) {
        return GRIB_BUFFER_TOO_SMALL;
    }

    memcpy(buffer,tmp,already_read);

    if((r->read(r->read_data,buffer+already_read,rest,&err) != rest) || err)
        return err;

    if(check7777 && !r->headers_only && (buffer[message_length-4] != '7' ||
            buffer[message_length-3] != '7' ||
            buffer[message_length-2] != '7' ||
            buffer[message_length-1] != '7')) {

        return GRIB_WRONG_LENGTH;
    }

    return GRIB_SUCCESS;
}

#define CHECK_TMP_SIZE(a) if(sizeof(tmp)<(a)) { fprintf(stderr,"%s:%d sizeof(tmp)<%s %d<%d\n", __FILE__,__LINE__,#a,(int)sizeof(tmp),(int)(a)); return GRIB_INTERNAL_ARRAY_TOO_SMALL; }

#define GROW_BUF_IF_REQUIRED(desired_length) if(buf->length<(desired_length)) { grib_grow_buffer(c, buf, desired_length);tmp=buf->data; }

#define UINT3(a,b,c) (size_t)((a<<16) + (b<<8) + c);

static int read_GRIB(reader* r)
{
    unsigned char *tmp=NULL;
    size_t length = 0;
    size_t total_length = 0;
    long edition = 0;
    int  err = 0;
    int i = 0 ,j;
    size_t sec1len = 0;
    size_t sec2len = 0;
    size_t sec3len = 0;
    size_t sec4len = 0;
    unsigned long flags;
    size_t buflen = 32768;   /* See ECC-515: was 16368 */
    grib_context* c;
    grib_buffer* buf;

    /*TODO proper context*/
    c=grib_context_get_default();
    tmp=(unsigned char*)malloc(buflen);
    if (!tmp)
        return GRIB_OUT_OF_MEMORY;
    buf=grib_new_buffer(c,tmp,buflen);
    buf->property = GRIB_MY_BUFFER;

    tmp[i++] = 'G';
    tmp[i++] = 'R';
    tmp[i++] = 'I';
    tmp[i++] = 'B';

    r->offset=r->tell(r->read_data)-4;

    if(r->read(r->read_data,&tmp[i],3,&err) != 3 || err)
        return err;

    length= UINT3(tmp[i],tmp[i+1],tmp[i+2]);
    i+=3;

    /* Edition number */
    if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
        return err;

    edition = tmp[i++];

    switch(edition)
    {
    case 1:
        if (r->headers_only) {
            /* Read section 1 length */
            if(r->read(r->read_data,&tmp[i],3,&err) != 3 || err)
                return err;

            sec1len=UINT3(tmp[i],tmp[i+1],tmp[i+2]);
            i+=3;
            /* Read section 1. 3 = length */
            if((r->read(r->read_data,tmp+i,sec1len-3,&err) != sec1len-3) || err)
                return err;
            flags = tmp[15];

            i += sec1len-3;

            GROW_BUF_IF_REQUIRED(i+3);

            if(flags & (1<<7)) {
                /* Section 2 */
                if(r->read(r->read_data,&tmp[i],3,&err) != 3 || err)
                    return err;

                sec2len=UINT3(tmp[i],tmp[i+1],tmp[i+2]);
                GROW_BUF_IF_REQUIRED(i+sec2len);
                i+=3;
                /* Read section 2 */
                if((r->read(r->read_data,tmp+i,sec2len-3,&err) != sec2len-3) || err)
                    return err;
                i += sec2len-3;
            }


            if(flags & (1<<6)) {

                /* Section 3 */
                GROW_BUF_IF_REQUIRED(i+3);
                for(j=0;j<3;j++)
                {
                    if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                        return err;

                    sec3len <<= 8;
                    sec3len |= tmp[i];
                    i++;
                }

                /* Read section 3 */
                GROW_BUF_IF_REQUIRED(i + sec3len);
                if((r->read(r->read_data,tmp+i,sec3len-3,&err) != sec3len-3) || err)
                    return err;
                i += sec3len-3;
            }

            GROW_BUF_IF_REQUIRED(i + 11);

            /* Section 4 */
            for(j=0;j<3;j++)
            {
                if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                    return err;

                sec4len <<= 8;
                sec4len |= tmp[i];
                i++;
            }

            /* we don't read the data, only headers */
            if((r->read(r->read_data,tmp+i,8,&err) != 8) || err)
                return err;

            i+=8;

            total_length=length;
            /* length=8+sec1len + sec2len+sec3len+11; */
            length=i;
            err=r->seek(r->read_data,total_length-length-1);

        }
        else if(length & 0x800000)
        {

            /* Large GRIBs */

            /* Read section 1 length */
            for(j=0;j<3;j++)
            {
                if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                    return err;

                sec1len <<= 8;
                sec1len |= tmp[i];
                i++;
            }

            /* table version */
            if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
            /* center */
            if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
            /* process */
            if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
            /* grid */
            if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
            /* flags */
            if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err) return err;
            flags = tmp[i++];

            /* fprintf(stderr," sec1len=%d i=%d flags=%x\n",sec1len,i,flags); */

            GROW_BUF_IF_REQUIRED(8+sec1len +  4 + 3);

            /* Read section 1. 3 = length, 5 = table,center,process,grid,flags */
            if((r->read(r->read_data,tmp+i,sec1len-3-5,&err) != sec1len-3-5) || err)
                return err;

            i += sec1len-3-5;

            if(flags & (1<<7)) {
                /* Section 2 */
                for(j=0;j<3;j++)
                {
                    if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                        return err;

                    sec2len <<= 8;
                    sec2len |= tmp[i];
                    i++;
                }
                /* Read section 2 */
                GROW_BUF_IF_REQUIRED(i+sec2len);
                if((r->read(r->read_data,tmp+i,sec2len-3,&err) != sec2len-3) || err)
                    return err;
                i += sec2len-3;
            }

            GROW_BUF_IF_REQUIRED(sec1len +  sec2len + 4 + 3);

            if(flags & (1<<6)) {

                /* Section 3 */
                for(j=0;j<3;j++)
                {
                    if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                        return err;

                    sec3len <<= 8;
                    sec3len |= tmp[i];
                    i++;
                }

                /* Read section 3 */
                GROW_BUF_IF_REQUIRED(sec1len + sec2len + sec3len + 4 + 3);
                if((r->read(r->read_data,tmp+i,sec3len-3,&err) != sec3len-3) || err)
                    return err;
                i += sec3len-3;
            }

            /* fprintf(stderr,"%s sec1len=%d i=%d\n",type,sec1len,i); */

            GROW_BUF_IF_REQUIRED(sec1len + sec2len + sec3len + 4 + 3);


            for(j=0;j<3;j++)
            {
                if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                    return err;

                sec4len <<= 8;
                sec4len |= tmp[i];
                i++;
            }

            if(sec4len < 120)
            {
                /* Special coding */
                length &= 0x7fffff;
                length *= 120;
                length -= sec4len;
                length += 4;
            }
            else
            {
                /* length is already set to the right value */
            }

        }
        break;

    case 2:
    case 3:
        length = 0;

        if(sizeof(long) >= 8) {
            for(j=0;j<8;j++)
            {
                if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                    return err;

                length <<= 8;
                length |= tmp[i];
                i++;
            }
        }
        else
        {
            /* Check if the length fits in a long */
            for(j=0;j<4;j++)
            {
                if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                    return err;

                length <<= 8;
                length |= tmp[i];
                i++;
            }

            if(length)
                return GRIB_MESSAGE_TOO_LARGE; /* Message too large */

            for(j=0;j<4;j++)
            {
                if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                    return err;

                length <<= 8;
                length |= tmp[i];
                i++;
            }
        }
        break;

    default:
        r->seek_from_start(r->read_data,r->offset+4);
        return GRIB_UNSUPPORTED_EDITION;
        break;
    }

    /* Assert(i <= buf->length); */
    err=read_the_rest(r, length, tmp, i, 1);
    if (err)
      r->seek_from_start(r->read_data,r->offset+4);

    grib_buffer_delete(c,buf);

    return err;
}

static int read_PSEUDO(reader *r,const char* type)
{
    unsigned char tmp[32]; /* Should be enough */
    size_t sec1len = 0;
    size_t sec4len = 0;
    int  err = 0;
    int i = 0, j = 0;

    Assert( strlen(type) == 4 );
    for(j = 0; j < 4; j++)
    {
        tmp[i] = type[i];
        i++;
    }

    for(j=0;j<3;j++)
    {
        if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
            return err;

        sec1len <<= 8;
        sec1len |= tmp[i];
        i++;
    }

    /* fprintf(stderr,"%s sec1len=%d i=%d\n",type,sec1len,i); */
    CHECK_TMP_SIZE(sec1len + 4 + 3 );

    /* Read sectoin1 */
    if((r->read(r->read_data,tmp+i,sec1len-3,&err) != sec1len-3) || err)
        return err;

    i += sec1len-3;

    for(j=0;j<3;j++)
    {
        if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
            return err;

        sec4len <<= 8;
        sec4len |= tmp[i];
        i++;
    }

    /* fprintf(stderr,"%s sec4len=%d i=%d l=%d\n",type,sec4len,i,4+sec1len+sec4len+4); */

    Assert(i <= sizeof(tmp));
    return read_the_rest(r,4+sec1len+sec4len+4,tmp,i, 1);
}


static int read_HDF5_offset(reader *r, int length, unsigned long* v, unsigned char *tmp, int* i) {
    unsigned char buf[8];
    int j, k;
    int err = 0;


    if( (r->read(r->read_data, buf, length, &err) != length) || err) {
        return err;
    }

    k = *i;
    for(j = 0; j < length; j++) {
        tmp[k++] = buf[j];
    }
    *i = k;

    *v = 0;
    for(j = length-1; j >= 0; j--) {
        *v <<= 8;
        *v |= buf[j];
    }

    return 0;
}

static int read_HDF5(reader *r)
{
    /* See: http://www.hdfgroup.org/HDF5/doc/H5.format.html#Superblock */
    unsigned char tmp[49]; /* Should be enough */
    unsigned char buf[4];

    unsigned char version_of_superblock,  size_of_offsets, size_of_lengths, consistency_flags;
    unsigned long base_address, superblock_extension_address, end_of_file_address;

    int i = 0, j;
    int err = 0;
    grib_context* c = grib_context_get_default();

    tmp[i++] = 137;
    tmp[i++] = 'H';
    tmp[i++] = 'D';
    tmp[i++] = 'F';

    if( (r->read(r->read_data, buf, 4, &err) != 4) || err) {
        return err;
    }

    if( !(buf[0] == '\r' && buf[1] == '\n' && buf[2] == 26 && buf[3] == '\n')) {
        /* Invalid magic, we should not use grib_context_log without a context */
        grib_context_log(c, GRIB_LOG_ERROR,"read_HDF5: invalid signature");
        return GRIB_INVALID_MESSAGE;
    }

    for(j = 0; j < 4; j++) {
        tmp[i++] = buf[j];
    }

    if( (r->read(r->read_data, &version_of_superblock, 1,  &err)  != 1) || err) {
        return err;
    }

    tmp[i++] = version_of_superblock;

    if (version_of_superblock == 2 || version_of_superblock == 3) {
        if( (r->read(r->read_data, &size_of_offsets, 1, &err) != 1) || err) {
            return err;
        }

        tmp[i++] = size_of_offsets;

        if(size_of_offsets > 8) {
            grib_context_log(c, GRIB_LOG_ERROR,"read_HDF5: invalid size_of_offsets: %ld, only <= 8 is supported", (long)size_of_offsets);
            return GRIB_NOT_IMPLEMENTED;
        }

        if( (r->read(r->read_data, &size_of_lengths, 1, &err) != 1) || err) {
            return err;
        }

        tmp[i++] = size_of_lengths;

        if( (r->read(r->read_data, &consistency_flags, 1, &err) != 1) || err) {
            return err;
        }

        tmp[i++] = consistency_flags;

        err = read_HDF5_offset(r, size_of_offsets, &base_address, tmp, &i);
        if(err) {
            return err;
        }

        err = read_HDF5_offset(r, size_of_offsets, &superblock_extension_address, tmp, &i);
        if(err) {
            return err;
        }

        err = read_HDF5_offset(r, size_of_offsets, &end_of_file_address, tmp, &i);
        if(err) {
            return err;
        }
    } else if (version_of_superblock == 0 || version_of_superblock == 1) {
        char skip[4];
        unsigned long file_free_space_info;
        unsigned char version_of_file_free_space, version_of_root_group_symbol_table, version_number_shared_header, ch;

        if( (r->read(r->read_data, &version_of_file_free_space, 1, &err) != 1) || err) return err;
        tmp[i++] = version_of_file_free_space;

        if( (r->read(r->read_data, &version_of_root_group_symbol_table, 1, &err) != 1) || err) return err;
        tmp[i++] = version_of_root_group_symbol_table;

        if( (r->read(r->read_data, &ch, 1, &err) != 1) || err) return err; /* reserved */
        tmp[i++] = ch;

        if( (r->read(r->read_data, &version_number_shared_header, 1, &err) != 1) || err) return err;
        tmp[i++] = version_number_shared_header;

        if( (r->read(r->read_data, &size_of_offsets, 1, &err) != 1) || err) return err;
        tmp[i++] = size_of_offsets;
        if (size_of_offsets > 8) {
            grib_context_log(c, GRIB_LOG_ERROR,"read_HDF5: invalid size_of_offsets: %ld, only <= 8 is supported", (long)size_of_offsets);
            return GRIB_NOT_IMPLEMENTED;
        }

        if( (r->read(r->read_data, &size_of_lengths, 1, &err) != 1) || err) return err;
        tmp[i++] = size_of_lengths;

        if( (r->read(r->read_data, &ch, 1, &err) != 1) || err) return err; /*reserved*/
        tmp[i++] = ch;

        if( (r->read(r->read_data, &skip, 4, &err) != 4) || err) return err; /* Group Leaf/Internal Node K: 4 bytes */
        tmp[i++] = skip[0]; tmp[i++] = skip[1]; tmp[i++] = skip[2]; tmp[i++] = skip[3];

        if( (r->read(r->read_data, &skip, 4, &err) != 4) || err) return err; /* consistency_flags: 4 bytes */
        tmp[i++] = skip[0]; tmp[i++] = skip[1]; tmp[i++] = skip[2]; tmp[i++] = skip[3];

        if (version_of_superblock == 1) {
            /* Indexed storage internal node K and reserved: only in version 1 of superblock */
            if( (r->read(r->read_data, &skip, 4, &err) != 4) || err) return err;
            tmp[i++] = skip[0]; tmp[i++] = skip[1]; tmp[i++] = skip[2]; tmp[i++] = skip[3];
        }

        err = read_HDF5_offset(r, size_of_offsets, &base_address, tmp, &i);
        if (err) return err;

        err = read_HDF5_offset(r, size_of_offsets, &file_free_space_info, tmp, &i);
        if (err) return err;

        err = read_HDF5_offset(r, size_of_offsets, &end_of_file_address, tmp, &i);
        if (err) return err;
    } else {
        grib_context_log(c, GRIB_LOG_ERROR,"read_HDF5: invalid version of superblock: %ld", (long)version_of_superblock);
        return GRIB_NOT_IMPLEMENTED;
    }

    Assert(i <= sizeof(tmp));
    return read_the_rest(r, end_of_file_address, tmp, i, 0);
}

static int read_WRAP(reader *r)
{
    /* See: http://www.hdfgroup.org/HDF5/doc/H5.format.html#Superblock */
    unsigned char tmp[36]; /* Should be enough */
    unsigned char buf[8];

    unsigned long long length = 0;

    int i = 0, j;
    int err = 0;

    tmp[i++] = 'W';
    tmp[i++] = 'R';
    tmp[i++] = 'A';
    tmp[i++] = 'P';

    if( (r->read(r->read_data, buf, 8, &err) != 8) || err) {
        printf("error\n");
        return err;
    }

    for(j = 0; j < 8; j++) {
        length <<= 8;
        length |= buf[j];
        tmp[i++] = buf[j];
    }

    Assert(i <= sizeof(tmp));
    return read_the_rest(r, length, tmp, i, 1);
}

static int read_BUFR(reader *r)
{
    /* unsigned char tmp[65536];*/ /* Should be enough */
    size_t length = 0;
    long edition = 0;
    int  err = 0;
    int i = 0 ,j;
    size_t buflen=2048;
    unsigned char *tmp=NULL;
    grib_context *c=NULL;
    grib_buffer* buf=NULL;

    /*TODO proper context*/
    c=grib_context_get_default();
    tmp=(unsigned char*)malloc(buflen);
    if (!tmp)
        return GRIB_OUT_OF_MEMORY;
    buf=grib_new_buffer(c,tmp,buflen);
    buf->property = GRIB_MY_BUFFER;
    r->offset=r->tell(r->read_data)-4;

    tmp[i++] = 'B';
    tmp[i++] = 'U';
    tmp[i++] = 'F';
    tmp[i++] = 'R';

    for(j=0;j<3;j++)
    {
        if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
            return err;

        length <<= 8;
        length |= tmp[i];
        i++;
    }

    if(length==0) {
        return GRIB_INVALID_MESSAGE;
    }

    /* Edition number */
    if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
        return err;

    edition = tmp[i++];

    /* Assert(edition != 1); */

    switch (edition) {
      case 0:
      case 1:
        {
        int n;
        size_t sec1len = 0;
        size_t sec2len = 0;
        size_t sec3len = 0;
        size_t sec4len = 0;
        unsigned long flags;

        sec1len = length;

        /* table version */
        if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
        /* center */
        if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
        /* update */
        if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
        /* flags */
        if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err) return err;
        flags = tmp[i++];


        GROW_BUF_IF_REQUIRED(sec1len +  4 + 3);

        /* Read section 1. 3 = length, 5 = table,center,process,flags */

        n = sec1len - 8; /* Just a guess */
        if((r->read(r->read_data,tmp+i,n,&err) != n) || err)
            return err;

        i += n;

        if(flags & (1<<7)) {
            /* Section 2 */
            for(j=0;j<3;j++)
            {
                if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                    return err;

                sec2len <<= 8;
                sec2len |= tmp[i];
                i++;
            }

            GROW_BUF_IF_REQUIRED(sec1len +  sec2len + 4 + 3);

            /* Read section 2 */
            if((r->read(r->read_data,tmp+i,sec2len-3,&err) != sec2len-3) || err)
                return err;
            i += sec2len-3;
        }


        /* Section 3 */
        for(j=0;j<3;j++)
        {
            if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                return err;

            sec3len <<= 8;
            sec3len |= tmp[i];
            i++;
        }

        GROW_BUF_IF_REQUIRED(sec1len +  sec2len + sec3len + 4 + 3);

        /* Read section 3 */
        if((r->read(r->read_data,tmp+i,sec3len-3,&err) != sec3len-3) || err)
            return err;
        i += sec3len-3;

        for(j=0;j<3;j++)
        {
            if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
                return err;

            sec4len <<= 8;
            sec4len |= tmp[i];
            i++;
        }

        /* fprintf(stderr," sec1len=%d sec2len=%d sec3len=%d sec4len=%d\n",sec1len, sec2len,sec3len,sec4len); */
        length = 4 + sec1len + sec2len + sec3len + sec4len + 4;
        /* fprintf(stderr,"length = %d i = %d\n",length,i); */
        }
        break;
      case 2:
      case 3:
      case 4:
        break;
      default :
        r->seek_from_start(r->read_data,r->offset+4);
        return GRIB_UNSUPPORTED_EDITION;
    }


    /* Assert(i <= sizeof(tmp)); */
    err=read_the_rest(r, length, tmp, i, 1);
    if (err)
      r->seek_from_start(r->read_data,r->offset+4);

    grib_buffer_delete(c,buf);

    return err;
}

static int _read_any(reader *r, int grib_ok, int bufr_ok, int hdf5_ok, int wrap_ok)
{
    unsigned char c;
    int err = 0;
    unsigned long magic = 0;

    while(r->read(r->read_data,&c,1,&err) == 1 && err == 0)
    {
        magic <<= 8;
        magic |= c;

        switch(magic & 0xffffffff)
        {
        case GRIB:
            if(grib_ok)
            {
                err =  read_GRIB(r);
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
            }
            break;

        case BUFR:
            if(bufr_ok)
            {
                err =  read_BUFR(r);
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
            }
            break;

        case HDF5:
            if(hdf5_ok)
            {
                err =  read_HDF5(r);
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
            }
            break;

        case WRAP:
            if(wrap_ok)
            {
                err =  read_WRAP(r);
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
            }
            break;

        case BUDG:
            if(grib_ok)
            {
                err =  read_PSEUDO(r,"BUDG");
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
            }
            break;
        case DIAG:
            if(grib_ok)
            {
                err =  read_PSEUDO(r,"DIAG");
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
            }
            break;
        case TIDE:
            if(grib_ok)
            {
                err =  read_PSEUDO(r,"TIDE");
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
            }
            break;
        }
    }

    return err;
}
static int read_any(reader *r, int grib_ok, int bufr_ok, int hdf5_ok, int wrap_ok)
{
    int result = 0;
    GRIB_MUTEX_INIT_ONCE(&once,&init);
    GRIB_MUTEX_LOCK(&mutex1);
    result = _read_any(r, grib_ok, bufr_ok, hdf5_ok, wrap_ok);
    GRIB_MUTEX_UNLOCK(&mutex1);
    return result;
}

static int read_any_gts(reader *r)
{
    unsigned char c;
    int err = 0;
    unsigned char* buffer=NULL;
    unsigned long magic = 0;
    unsigned long start = 0x010d0d0a;
    unsigned long theEnd = 0x0d0d0a03;
    unsigned char tmp[128]={0,}; /* See ECC-735 */
    size_t message_size=0;
    size_t already_read=0;
    int i=0;

    while(r->read(r->read_data,&c,1,&err) == 1 && err == 0)
    {
        magic <<= 8;
        magic |= c;
        magic &= 0xffffffff;

        if (magic == start) {
            tmp[i++]=0x01;
            tmp[i++]=0x0d;
            tmp[i++]=0x0d;
            tmp[i++]=0x0a;

            r->offset=r->tell(r->read_data)-4;

            if(r->read(r->read_data,&tmp[i],6,&err) != 6 || err)
                return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */

            if (tmp[7] != 0x0d || tmp[8]!= 0x0d || tmp[9]!=0x0a) {
                r->seek(r->read_data,-6);
                continue;
            }
            magic=0;
            already_read=10;
            message_size=already_read;
            while(r->read(r->read_data,&c,1,&err) == 1 && err == 0) {
                message_size++;
                magic <<= 8;
                magic |= c;
                magic &= 0xffffffff;
                if (magic == theEnd) {
                    r->seek(r->read_data,already_read-message_size);
                    buffer = (unsigned char*)r->alloc(r->alloc_data,&message_size,&err);
                    if (!buffer) return GRIB_OUT_OF_MEMORY;
                    if (err) return err;
                    memcpy(buffer,tmp,already_read);
                    r->read(r->read_data,buffer+already_read,message_size-already_read,&err);
                    r->message_size=message_size;
                    return err;
                }
            }
        }
    }

    return err;
}

static int read_any_taf(reader *r)
{
    unsigned char c;
    int err = 0;
    unsigned char* buffer=NULL;
    unsigned long magic = 0;
    unsigned long start = 0x54414620;
    unsigned char tmp[1000]={0,}; /* Should be enough */
    size_t message_size=0;
    size_t already_read=0;
    int i=0;

    while(r->read(r->read_data,&c,1,&err) == 1 && err == 0)
    {
        magic <<= 8;
        magic |= c;
        magic &= 0xffffffff;

        if (magic == start) {
            tmp[i++]=0x54;
            tmp[i++]=0x41;
            tmp[i++]=0x46;
            tmp[i++]=0x20;

            r->offset=r->tell(r->read_data)-4;

            already_read=4;
            message_size=already_read;
            while(r->read(r->read_data,&c,1,&err) == 1 && err == 0) {
                message_size++;
                if (c == '=') {
                    r->seek(r->read_data,already_read-message_size);
                    buffer = (unsigned char*)r->alloc(r->alloc_data,&message_size,&err);
                    if (!buffer) return GRIB_OUT_OF_MEMORY;
                    if (err) return err;
                    memcpy(buffer,tmp,already_read);
                    r->read(r->read_data,buffer+already_read,message_size-already_read,&err);
                    r->message_size=message_size;
                    return err;
                }
            }
        }
    }

    return err;
}

static int read_any_metar(reader *r)
{
    unsigned char c;
    int err = 0;
    unsigned char* buffer=NULL;
    unsigned long magic = 0;
    unsigned long start = 0x4d455441;
    unsigned char tmp[32]={0,}; /* Should be enough */
    size_t message_size=0;
    size_t already_read=0;
    int i=0;

    while(r->read(r->read_data,&c,1,&err) == 1 && err == 0)
    {
        magic <<= 8;
        magic |= c;
        magic &= 0xffffffff;

        if (magic == start) {
            if (r->read(r->read_data,&c,1,&err) != 1 || err!=0)
                break;
            if (c == 'R' ) {
                tmp[i++]=0x4d;
                tmp[i++]=0x45;
                tmp[i++]=0x54;
                tmp[i++]=0x41;
                tmp[i++]='R';

                r->offset=r->tell(r->read_data)-4;

                already_read=5;
                message_size=already_read;
                while(r->read(r->read_data,&c,1,&err) == 1 && err == 0) {
                    message_size++;
                    if (c == '=') {
                        r->seek(r->read_data,already_read-message_size);
                        buffer = (unsigned char*)r->alloc(r->alloc_data,&message_size,&err);
                        if (!buffer) return GRIB_OUT_OF_MEMORY;
                        if (err) return err;
                        memcpy(buffer,tmp,already_read);
                        r->read(r->read_data,buffer+already_read,message_size-already_read,&err);
                        r->message_size=message_size;
                        return err;
                    }
                }
            }
        }
    }

    return err;
}

off_t stdio_tell(void* data)
{
    FILE* f = (FILE*)data;
    return ftello(f);
}

int stdio_seek(void* data,off_t len)
{
    FILE* f = (FILE*)data;
    int err=0;
    if (fseeko(f,len,SEEK_CUR)) err=GRIB_IO_PROBLEM;
    return err;
}

int stdio_seek_from_start(void* data,off_t len)
{
    FILE* f = (FILE*)data;
    int err=0;
    if (fseeko(f,len,SEEK_SET)) err=GRIB_IO_PROBLEM;
    return err;
}

size_t stdio_read(void* data,void* buf,size_t len,int* err)
{
    FILE* f = (FILE*)data;
    size_t n;
    /* char iobuf[1024*1024]; */

    if (len==0) return 0;

    /* setvbuf(f,iobuf,_IOFBF,sizeof(iobuf)); */
    n   = fread(buf,1,len,f);
    /* fprintf(stderr,"read %d = %x %c\n",1,(int)buf[0],buf[0]); */
    if(n != len) {
        /* fprintf(stderr,"Failed to read %d, only got %d\n",len,n); */
        *err               = GRIB_IO_PROBLEM;
        if(feof(f))   *err = GRIB_END_OF_FILE;
        if(ferror(f)) *err = GRIB_IO_PROBLEM;
    }
    return n;
}

/*================== */
typedef struct user_buffer {
    void*    user_buffer;
    size_t   buffer_size;
} user_buffer;

static void* user_provider_buffer(void *data,size_t* length,int *err)
{
    user_buffer *u  = (user_buffer*)data;
    *length = u->buffer_size;
    return u->user_buffer;
}

static
int _wmo_read_any_from_file(FILE* f,void* buffer,size_t* len,int grib_ok,int bufr_ok, int hdf5_ok, int wrap_ok)
{
    int         err;
    user_buffer u;
    reader      r;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.message_size    = 0;
    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 0;

    err            = read_any(&r, grib_ok, bufr_ok, hdf5_ok, wrap_ok);
    *len           = r.message_size;

    return err;
}

int wmo_read_any_from_file(FILE* f,void* buffer,size_t* len)
{
    return _wmo_read_any_from_file(f, buffer, len, 1, 1, 1, 1);
}

int wmo_read_grib_from_file(FILE* f,void* buffer,size_t* len)
{
    return _wmo_read_any_from_file(f, buffer, len, 1, 0, 0, 0);
}

int wmo_read_bufr_from_file(FILE* f,void* buffer,size_t* len)
{
    return _wmo_read_any_from_file(f, buffer, len, 0, 1, 0, 0);
}

int wmo_read_gts_from_file(FILE* f,void* buffer,size_t* len)
{
    int         err;
    user_buffer u;
    reader      r;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.message_size    = 0;
    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 0;

    err            = read_any_gts(&r);
    *len           = r.message_size;

    return err;
}

int wmo_read_taf_from_file(FILE* f,void* buffer,size_t* len)
{
    int         err;
    user_buffer u;
    reader      r;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 0;

    err            = read_any_taf(&r);
    *len           = r.message_size;

    return err;
}

int wmo_read_metar_from_file(FILE* f,void* buffer,size_t* len)
{
    int         err;
    user_buffer u;
    reader      r;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 0;

    err            = read_any_metar(&r);
    *len           = r.message_size;

    return err;
}

/*================== */

typedef struct stream_struct {

    void* stream_data;
    long (*stream_proc)(void*,void* buffer,long len);

} stream_struct;

static off_t stream_tell(void* data)
{
    return 0;
}

static int stream_seek(void* data,off_t len)
{
    return 0;
}
static size_t stream_read(void* data,void* buffer,size_t len,int* err)
{
    stream_struct *s = (stream_struct*)data;
    long n = len;

    if(n != len) {
        /* size_t cannot be coded into long */
        *err  = GRIB_INTERNAL_ERROR;
        return -1;
    }

    n =  s->stream_proc(s->stream_data,buffer,len);
    if(n != len) {
        *err  = GRIB_IO_PROBLEM;
        if(n == -1) *err = GRIB_END_OF_FILE;
    }
    return n;
}

/*================== */


static void* allocate_buffer(void *data,size_t* length,int *err)
{
    alloc_buffer *u  = (alloc_buffer*)data;
    u->buffer = malloc(*length);
    u->size=*length;
    if(u->buffer == NULL)
        *err = GRIB_OUT_OF_MEMORY; /* Cannot allocate buffer */
    return u->buffer;
}

int wmo_read_any_from_stream(void* stream_data,long (*stream_proc)(void*,void* buffer,long len) ,void* buffer,size_t* len)
{
    int           err;
    stream_struct s;
    user_buffer   u;
    reader        r;

    s.stream_data = stream_data;
    s.stream_proc = stream_proc;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.message_size    = 0;
    r.offset          = 0;
    r.read_data       = &s;
    r.read            = &stream_read;
    r.seek            = &stream_seek;
    r.seek_from_start = &stream_seek;
    r.tell            = &stream_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 0;

    err            = read_any(&r, 1, 1, 1, 1);
    *len           = r.message_size;

    return err;
}

void* wmo_read_any_from_stream_malloc(void* stream_data,long (*stream_proc)(void*,void* buffer,long len) ,size_t *size, int* err)
{
    alloc_buffer  u;
    stream_struct s;
    reader        r;

    u.buffer = NULL;

    s.stream_data = stream_data;
    s.stream_proc = stream_proc;

    r.message_size    = 0;
    r.offset          = 0;
    r.read_data       = &s;
    r.read            = &stream_read;
    r.seek            = &stream_seek;
    r.seek_from_start = &stream_seek;
    r.tell            = &stream_tell;
    r.alloc_data      = &u;
    r.alloc           = &allocate_buffer;
    r.headers_only    = 0;

    *err           = read_any(&r, 1, 1, 1, 1);
    *size          = r.message_size;

    return u.buffer;
}

/*================== */


void *wmo_read_gts_from_file_malloc(FILE* f,int headers_only,size_t *size,off_t *offset,int* err)
{
    alloc_buffer u;
    reader       r;

    u.buffer       = NULL;
    r.offset       = 0;

    r.message_size    = 0;
    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &allocate_buffer;
    r.headers_only    = headers_only;

    *err           = read_any_gts(&r);
    *size          = r.message_size;
    *offset        = r.offset;

    return u.buffer;
}

void *wmo_read_taf_from_file_malloc(FILE* f,int headers_only,size_t *size,off_t *offset,int* err)
{
    alloc_buffer u;
    reader       r;

    u.buffer       = NULL;

    r.offset          = 0;
    r.message_size    = 0;
    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &allocate_buffer;
    r.headers_only    = headers_only;

    *err           = read_any_taf(&r);
    *size          = r.message_size;
    *offset        = r.offset;

    return u.buffer;
}

void *wmo_read_metar_from_file_malloc(FILE* f,int headers_only,size_t *size,off_t *offset,int* err)
{
    alloc_buffer u;
    reader       r;

    u.buffer       = NULL;

    r.message_size    = 0;
    r.read_data       = f;
    r.offset          = 0;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &allocate_buffer;
    r.headers_only    = headers_only;

    *err           = read_any_metar(&r);
    *size          = r.message_size;
    *offset        = r.offset;

    return u.buffer;
}

static void *_wmo_read_any_from_file_malloc(FILE* f,int* err,size_t *size,off_t *offset,
        int grib_ok,int bufr_ok, int hdf5_ok, int wrap_ok, int headers_only)
{
    alloc_buffer u;
    reader       r;

    u.buffer       = NULL;
    u.size         = 0;

    r.message_size    = 0;
    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &allocate_buffer;
    r.headers_only    = headers_only;
    r.offset          = 0;

    *err           = read_any(&r, grib_ok, bufr_ok, hdf5_ok, wrap_ok);

    *size          = r.message_size;
    *offset        = r.offset;

    return u.buffer;
}
void *wmo_read_any_from_file_malloc(FILE* f,int headers_only,size_t *size,off_t *offset,int* err)
{
    return _wmo_read_any_from_file_malloc(f,err,size,offset, 1, 1, 1, 1, headers_only);
}

void *wmo_read_grib_from_file_malloc(FILE* f,int headers_only,size_t *size,off_t *offset,int* err)
{
    return _wmo_read_any_from_file_malloc(f,err,size,offset, 1, 0, 0, 0, headers_only);
}

void *wmo_read_bufr_from_file_malloc(FILE* f,int headers_only,size_t *size,off_t *offset,int* err)
{
    return _wmo_read_any_from_file_malloc(f,err,size,offset, 0, 1, 0, 0, headers_only);
}


/* ======================================= */

typedef struct context_alloc_buffer {
    grib_context*  ctx;
    void*          buffer;
    size_t         length;
} context_alloc_buffer;

static void* context_allocate_buffer(void *data,size_t* length,int *err)
{
    context_alloc_buffer *u  = (context_alloc_buffer*)data;
    u->buffer = grib_context_malloc(u->ctx,*length);
    u->length = *length;

    if(u->buffer == NULL)
        *err = GRIB_OUT_OF_MEMORY; /* Cannot allocate buffer */
    return u->buffer;
}


int grib_read_any_headers_only_from_file(grib_context* ctx,FILE* f,void* buffer,size_t* len)
{
    int         err;
    user_buffer u;
    reader      r;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.message_size    = 0;
    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 1;

    err            = read_any(&r, 1, GRIB_API_READS_BUFR, GRIB_API_READS_HDF5, GRIB_API_READS_WRAP);

    *len           = r.message_size;

    return err;
}

int grib_read_any_from_file(grib_context* ctx,FILE* f,void* buffer,size_t* len)
{
    int         err;
    user_buffer u;
    reader      r;
    off_t       offset;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.message_size    = 0;
    r.read_data       = f;
    r.read            = &stdio_read;
    r.seek            = &stdio_seek;
    r.seek_from_start = &stdio_seek_from_start;
    r.tell            = &stdio_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 0;

    offset=ftello(f);

    err            = read_any(&r, 1, GRIB_API_READS_BUFR, GRIB_API_READS_HDF5, GRIB_API_READS_WRAP);

    if (err==GRIB_BUFFER_TOO_SMALL) {
        if (fseeko(f,offset,SEEK_SET))
            err=GRIB_IO_PROBLEM;
    }

    *len           = r.message_size;

    return err;
}

/* ======================================= */

typedef struct memory_read_data {
    unsigned char  *data;
    size_t          data_len;
} memory_read_data;

static off_t memory_tell(void* data)
{
    return 0;
}

static int memory_seek(void* data,off_t len)
{
    return 0;
}

static size_t memory_read(void* data,void* buf,size_t len,int* err)
{
    memory_read_data *m = (memory_read_data*)data;

    if(len == 0)
    {
        *err = GRIB_END_OF_FILE;
        return 0;
    }
    else {
        size_t l = len > m->data_len ? m->data_len : len;
        memcpy(buf,m->data,l);
        m->data_len -= l;
        m->data     += l;
        return l;
    }
}

int grib_read_any_from_memory_alloc(grib_context* ctx,unsigned char** data,size_t* data_length,void **buffer,size_t* length)
{
    int err;
    memory_read_data     m;
    context_alloc_buffer u;
    reader       r;

    m.data         = *data;
    m.data_len     = *data_length;

    u.buffer       = NULL;
    u.length       = 0;
    u.ctx          = ctx ? ctx : grib_context_get_default();

    r.message_size    = 0;
    r.read_data       = &m;
    r.read            = &memory_read;
    r.seek            = &memory_seek;
    r.seek_from_start = &memory_seek;
    r.tell            = &memory_tell;
    r.alloc_data      = &u;
    r.alloc           = &context_allocate_buffer;
    r.headers_only    = 0;

    err            = read_any(&r, 1, GRIB_API_READS_BUFR, GRIB_API_READS_HDF5, GRIB_API_READS_WRAP);
    *buffer        = u.buffer;
    *length        = u.length;

    *data_length   = m.data_len;
    *data          = m.data;

    return err;
}

int grib_read_any_from_memory(grib_context* ctx,unsigned char** data,size_t* data_length,void* buffer,size_t* len)
{
    int         err;
    memory_read_data     m;
    user_buffer u;
    reader      r;

    m.data         = *data;
    m.data_len     = *data_length;

    u.user_buffer  = buffer;
    u.buffer_size  = *len;

    r.message_size    = 0;
    r.read_data       = &m;
    r.read            = &memory_read;
    r.seek            = &memory_seek;
    r.seek_from_start = &memory_seek;
    r.tell            = &memory_tell;
    r.alloc_data      = &u;
    r.alloc           = &user_provider_buffer;
    r.headers_only    = 0;

    err            = read_any(&r, 1, GRIB_API_READS_BUFR, GRIB_API_READS_HDF5, GRIB_API_READS_WRAP);
    *len           = r.message_size;

    *data_length   = m.data_len;
    *data          = m.data;

    return err;
}

int grib_count_in_file(grib_context* c, FILE* f,int* n)
{
    int err=0;
    *n=0;
    if (!c) c=grib_context_get_default();

    if ( c->multi_support_on )
    {
        /* GRIB-395 */
        grib_handle* h=NULL;
        while ((h=grib_handle_new_from_file(c, f , &err))!=NULL) {
            grib_handle_delete(h);
            (*n)++;
        }
    }
    else
    {
        void* mesg=NULL;
        size_t size=0;
        off_t offset=0;
        while ( (mesg=wmo_read_any_from_file_malloc ( f,0,  &size,&offset,&err))!=NULL && err==GRIB_SUCCESS) {
            grib_context_free(c,mesg);
            (*n)++;
        }
    }

    rewind(f);

    return err==GRIB_END_OF_FILE ? 0 : err;
}

int grib_count_in_filename(grib_context* c, const char* filename, int* n)
{
    int err=0;
    FILE* fp = NULL;
    if (!c) c=grib_context_get_default();
    fp = fopen(filename, "r");
    if (!fp) {
        grib_context_log(c, GRIB_LOG_ERROR,"grib_count_in_filename: Unable to read file \"%s\"", filename);
        perror(filename);
        return GRIB_IO_PROBLEM;
    }
    err = grib_count_in_file(c, fp, n);
    fclose(fp);
    return err;
}

/*****************************************************************************
 *
 * This following code implements a buffered I/O interface to URL reads using
 * libcurl. Based on example "fopen.c" in libcurl sources, with additions to
 * support fseek using range requests.
 *
 * Copyright (c) 2003, 2017 Simtec Electronics
 *
 * Re-implemented by Vincent Sanders <vince@kyllikki.org> with extensive
 * reference to original curl example code
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * This example requires libcurl 7.9.7 or later.
 */

#ifdef HAVE_CURL

#ifdef _WIN32
#define WAITMS(x) Sleep(x)
#else
/* Portable sleep for platforms other than Windows. */
#define WAITMS(x)                               \
  struct timeval wait = { 0, (x) * 1000 };      \
  (void)select(0, NULL, NULL, NULL, &wait);
#endif // _WIN32
 
typedef struct url_struct {
    CURL *handle;               /* libcurl handle */
    char *buffer;               /* buffer to store cached data */
    size_t buffer_len;          /* currently allocated buffers length */
    size_t buffer_pos;          /* end of data in buffer */
    off_t stream_len;           /* content length */
    off_t stream_pos;           /* current position in content */
    int still_running;          /* background url fetch still in progress */
    int allow_range;            /* server allows range transfers */
} url_struct;

static CURLM *url_multi_handle;

/* called when header data is received */
static size_t url_header_callback(char *buffer, size_t size, size_t nitems, void *userp)
{
    url_struct *us = (url_struct *)userp;
    size_t realsize = size * nitems;

    /* determine if range transfers are allowed */
    const char *accept_line = "Accept-Ranges: bytes";
    if (realsize >= strlen(accept_line)
        && strncmp(buffer, accept_line, strlen(accept_line)) == 0) {
        us->allow_range = 1;
    }
  
    return realsize;
}

/* called when data is received */
static size_t url_write_callback(char *buffer, size_t size, size_t nitems, void *userp)
{
  char *newbuff;
  size_t rembuff;

  url_struct *us = (url_struct *)userp;
  size_t realsize = size * nitems;

  rembuff = us->buffer_len - us->buffer_pos; /* remaining space in buffer */

  if (realsize > rembuff) {
    /* not enough space in buffer */
    newbuff = realloc(us->buffer, us->buffer_len + (realsize - rembuff));
    if (newbuff == NULL) {
      /* fprintf(stderr, "callback buffer grow failed\n"); */
      realsize = rembuff;
    }
    else {
      /* realloc succeeded - increase buffer size*/
      us->buffer_len += realsize - rembuff;
      us->buffer = newbuff;
    }
  }

  memcpy(&us->buffer[us->buffer_pos], buffer, realsize);
  us->buffer_pos += realsize;
  
  return realsize;
}

/* attempt to fill the read buffer up to requested number of bytes */
static int url_fill_buffer(url_struct *us, size_t n)
{
    fd_set fdread;
    fd_set fdwrite;
    fd_set fdexcep;
    struct timeval timeout;
    int rc;
    CURLMcode mc; /* curl_multi_fdset() return code */

    /* only attempt to fill buffer if transactions still running and buffer
     * doesn't exceed required size already
     */
    if ((!us->still_running) || (us->buffer_pos > n))
        return 0;

    /* attempt to fill buffer */
    do {
        int maxfd = -1;
        long curl_timeo = -1;

        FD_ZERO(&fdread);
        FD_ZERO(&fdwrite);
        FD_ZERO(&fdexcep);

        /* set a suitable timeout to fail on */
        timeout.tv_sec = 30; /* 30 seconds */
        timeout.tv_usec = 0;

        curl_multi_timeout(url_multi_handle, &curl_timeo);
        if (curl_timeo >= 0) {
            timeout.tv_sec = curl_timeo / 1000;
            if (timeout.tv_sec > 1)
                timeout.tv_sec = 1;
            else
                timeout.tv_usec = (curl_timeo % 1000) * 1000;
        }

        /* get file descriptors from the transfers */
        mc = curl_multi_fdset(url_multi_handle, &fdread, &fdwrite, &fdexcep, &maxfd);

        if (mc != CURLM_OK) {
            /* fprintf(stderr, "curl_multi_fdset() failed, code %d.\n", mc); */
            break;
        }

        /* On success the value of maxfd is guaranteed to be >= -1. We call
         * select(maxfd + 1, ...); specially in case of (maxfd == -1) there are
         * no fds ready yet so we call select(0, ...) --or Sleep() on Windows--
         * to sleep 100ms, which is the minimum suggested value in the
         * curl_multi_fdset() doc.
         */
        if (maxfd == -1) {
            WAITMS(100);
        }
        else {
            /* Note that on some platforms 'timeout' may be modified by select().
               If you need access to the original value save a copy beforehand. */
            rc = select(maxfd + 1, &fdread, &fdwrite, &fdexcep, &timeout);
        }

        switch(rc) {
        case -1:
            /* select error */
            break;
        case 0:
        default:
            /* timeout or readable/writable sockets */
            curl_multi_perform(url_multi_handle, &us->still_running);
            break;
        }
    } while (us->still_running && (us->buffer_pos < n));
    
    return 1;
}

/* remove n bytes from the front of a buffer */
static int url_use_buffer(url_struct *us, size_t n)
{
    if ((us->buffer_pos - n) <= 0) {
        /* clear buffer - write will recreate */
        free(us->buffer);
        us->buffer = NULL;
        us->buffer_pos = 0;
        us->buffer_len = 0;
    }
    else {
        /* move remaining bytes to the front */
        memmove(us->buffer,
                &us->buffer[n],
                (us->buffer_pos - n));

        us->buffer_pos -= n;
    }
    return 0;
}

url_struct *url_open(const char *url)
{
    url_struct *us;

    us = calloc(1, sizeof(url_struct));
    if (!us)
        return NULL;

    us->handle = curl_easy_init();
    
    /* send HEAD request to initialize header data */
    curl_easy_setopt(us->handle, CURLOPT_URL, url);
    curl_easy_setopt(us->handle, CURLOPT_VERBOSE, 0L);
    curl_easy_setopt(us->handle, CURLOPT_NOBODY, 1L);
    curl_easy_setopt(us->handle, CURLOPT_HEADERDATA, us);
    curl_easy_setopt(us->handle, CURLOPT_HEADERFUNCTION, url_header_callback);
    
    CURLcode rc;
    rc = curl_easy_perform(us->handle);
    if (rc == CURLE_OK) {
        curl_off_t cl;
        curl_easy_getinfo(us->handle, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &cl);
        if (cl > 0) {
            us->stream_len = (off_t)cl;
        }
    }
    
    /* reset the handle for data transfer */
    curl_easy_reset(us->handle);
    curl_easy_setopt(us->handle, CURLOPT_URL, url);
    curl_easy_setopt(us->handle, CURLOPT_VERBOSE, 0L);
    curl_easy_setopt(us->handle, CURLOPT_WRITEDATA, us);
    curl_easy_setopt(us->handle, CURLOPT_WRITEFUNCTION, url_write_callback);
    
    /* initialize the global multi handle on first call */
    if (!url_multi_handle)
        url_multi_handle = curl_multi_init();
    
    curl_multi_add_handle(url_multi_handle, us->handle);

    /* start non-blocking transfers */
    curl_multi_perform(url_multi_handle, &us->still_running);
    
    /* if still_running is 0 now, cleanup and set NULL */
    if ((us->buffer_pos == 0) && (!us->still_running)) {
        curl_multi_remove_handle(url_multi_handle, us->handle);
        curl_easy_cleanup(us->handle);
        free(us);
        us = NULL;
    }
  
    return us;
}

int url_close(url_struct *us)
{
    /* wait for pending transfers */
    url_fill_buffer(us, 0L);
    if (us->still_running)
        return 1;

    /* make sure the easy handle is not in the multi handle anymore */
    curl_multi_remove_handle(url_multi_handle, us->handle);

    /* cleanup */
    curl_easy_cleanup(us->handle);
    free(us->buffer);
    free(us);

    return 0; /* OK */
}

void url_rewind(url_struct *us)
{
    /* halt transaction */
    curl_multi_remove_handle(url_multi_handle, us->handle);

    /* clear buffer and stream pos */
    free(us->buffer);
    us->buffer = NULL;
    us->buffer_pos = 0;
    us->buffer_len = 0;
    us->stream_pos = 0;

    /* set range to beginning */
    if (us->allow_range) {
        curl_easy_setopt(us->handle, CURLOPT_RESUME_FROM_LARGE, 0L); 
    }

    /* restart transaction */
    curl_multi_add_handle(url_multi_handle, us->handle);
}

int url_eof(url_struct *us)
{
    if ((us->buffer_pos == 0) && (!us->still_running))
        return 1;
    
    return 0; /* OK */
}

off_t url_tell(void* data)
{
    url_struct *us = (url_struct*)data;
    if (us == NULL) {
        errno = GRIB_IO_PROBLEM;
        return -1;
    }
    return us->stream_pos;
}

int url_seek_internal(url_struct *us, off_t offset, int whence)
{
    /* calculate the new start position */
    curl_off_t pos = us->stream_pos;

    if (whence == SEEK_SET) {
        pos = offset;
    }
    else if (whence == SEEK_CUR) {
        pos = us->stream_pos + offset;
    }
    else if (whence == SEEK_END) {
        pos = us->stream_len + offset;
    }
    
    /* ensure position is positive */
    if (pos < 0) {
        errno = EINVAL;
        return -1;
    }
    
    if (us->allow_range) {
        /* halt transaction */
        curl_multi_remove_handle(url_multi_handle, us->handle);

        /* clear the buffer */
        free(us->buffer);
        us->buffer = NULL;
        us->buffer_pos = 0;
        us->buffer_len = 0;

        /* enable range transfer */
        curl_easy_setopt(us->handle, CURLOPT_RESUME_FROM_LARGE, pos);

        /* restart transaction */
        curl_multi_add_handle(url_multi_handle, us->handle);
    }
    else {
        /* range transfer disabled */
        /* TODO: reuse existing buffer to improve efficiency */
        url_rewind(us);
        url_fill_buffer(us, (size_t)pos);
        url_use_buffer(us, (size_t)pos);  
    }
    
    us->stream_pos = (off_t)pos;

    return 0; /* OK */
}

int url_seek(void *data, off_t len)
{
    url_struct *us = (url_struct*)data;
    int err = 0;
    if (url_seek_internal(us, len, SEEK_CUR)) err = GRIB_IO_PROBLEM;
    return err;
}

int url_seek_from_start(void *data, off_t len)
{
    url_struct *us = (url_struct*)data;
    int err=0;
    if (url_seek_internal(us, len, SEEK_SET)) err = GRIB_IO_PROBLEM;
    return err;
}

size_t url_read_internal(void *ptr, size_t size, size_t count, url_struct *us)
{
    size_t realsize = count * size;
    
    url_fill_buffer(us, realsize);

    /* check if there's data in the buffer - if not fill_buffer()
     * either errored or EOF */
    if (!us->buffer_pos)
    return 0;

    /* ensure only available data is considered */
    if (us->buffer_pos < realsize)
        realsize = us->buffer_pos;

    /* xfer data to caller */
    if (ptr)
        memcpy(ptr, us->buffer, realsize);

    url_use_buffer(us, realsize);
    us->stream_pos += realsize;

    size_t nitems = realsize / size; /* number of items */

    return nitems;
}

size_t url_read(void *data, void *buf, size_t len, int *err)
{
    url_struct *us = (url_struct*)data;
    size_t n;

    if (len == 0) return 0;

    n = url_read_internal(buf, 1, len, us);
    
    if (n != len) {
        if (url_eof(us))
            *err = GRIB_END_OF_FILE;
        else
            *err = GRIB_IO_PROBLEM;
    }
    
    return n;
}

/*================== */

static void *_wmo_read_any_from_url_malloc(void *us, int *err, size_t *size, off_t *offset,
        int grib_ok, int bufr_ok, int hdf5_ok, int wrap_ok, int headers_only)
{
    alloc_buffer u;
    reader       r;

    u.buffer       = NULL;
    u.size         = 0;

    r.message_size    = 0;
    r.read_data       = us;
    r.read            = &url_read;
    r.seek            = &url_seek;
    r.seek_from_start = &url_seek_from_start;
    r.tell            = &url_tell;
    r.alloc_data      = &u;
    r.alloc           = &allocate_buffer;
    r.headers_only    = headers_only;
    r.offset          = 0;

    *err           = read_any(&r, grib_ok, bufr_ok, hdf5_ok, wrap_ok);

    *size          = r.message_size;
    *offset        = r.offset;

    return u.buffer;
}

void *wmo_read_any_from_url_malloc(void *us, int headers_only, size_t *size, off_t *offset, int* err)
{
    return _wmo_read_any_from_url_malloc(us, err, size, offset, 1, 1, 1, 1, headers_only);
}

void *wmo_read_grib_from_url_malloc(void *us, int headers_only, size_t *size, off_t *offset, int* err)
{
    return _wmo_read_any_from_url_malloc(us, err, size, offset, 1, 0, 0, 0, headers_only);
}

void *wmo_read_bufr_from_url_malloc(void *us, int headers_only, size_t *size, off_t *offset, int* err)
{
    return _wmo_read_any_from_url_malloc(us, err, size, offset, 0, 1, 0, 0, headers_only);
}

#endif // HAVE_CURL
