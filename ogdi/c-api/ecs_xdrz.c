#include <rpc/types.h>
#include <rpc/xdr.h>

#include "ecs.h"
#include "zlib.h"

/* Variables used for compression/decompression routines */
static void *obuf = NULL;
static unsigned int obufsize = 0;
static int enough = 0;

extern bool_t xdr_ecs_Result_Work(XDR *xdrs, ecs_Result *objp);

bool_t
xdr_ecs_Result_Free(XDR *xdrs, ecs_Result *objp)
{
	if (objp->compression.ctype == ECS_COMPRESS_NONE) {
		return xdr_ecs_Result_Work(xdrs, objp);
	}
	if (objp->compression.ctype != ECS_COMPRESS_ZLIB) {
		return (FALSE);
	}
	if (objp->compression.cblksize != 0) {
		return xdr_ecs_Result_Work(xdrs, objp);
	}
	return (TRUE);
}

bool_t
xdr_ecs_Result_Decode(XDR *xdrs, ecs_Result *objp)
{
	XDR mem_xdrs;
	void *zbuf;
	z_stream z;
	int status;
	int size;

	if (!xdr_u_int(xdrs, &objp->compression.cfullsize)) {
		return FALSE;
	}
	if (objp->compression.ctype == ECS_COMPRESS_NONE) {
		return xdr_ecs_Result_Work(xdrs, objp);
	}
	if (objp->compression.ctype != ECS_COMPRESS_ZLIB) {
		return (FALSE);
	}
	if (objp->compression.cblksize == 0) {
		return xdr_ecs_Result_Work(xdrs, objp);
	}

	if (obufsize < objp->compression.cfullsize) {
		free(obuf);
		obuf = malloc(objp->compression.cfullsize);
		if (obuf == NULL) {
			obufsize = 0;
			return (FALSE);
		}
		obufsize = objp->compression.cfullsize;
	}

	zbuf = malloc(objp->compression.cblksize);
	if (zbuf == NULL) {
		return (FALSE);
	}

	/*
	 * First, read in the compressed data and decompress into obuf
	 */
	xdrmem_create(&mem_xdrs, obuf, objp->compression.cfullsize, XDR_DECODE);
	z.zalloc = NULL;
	z.zfree = NULL;
	z.opaque = NULL;
	if (inflateInit(&z) != Z_OK) {
		free(zbuf);
		return (FALSE);
	}
	z.next_out = obuf;
	z.avail_out = objp->compression.cfullsize;
	do {
		if (! xdr_bytes(xdrs, (char **) &zbuf, &size,
				objp->compression.cblksize))
		{
			xdr_destroy(&mem_xdrs);
			return (FALSE);
		}
		z.next_in = zbuf;
		z.avail_in = objp->compression.cblksize;
		status = inflate(&z, Z_NO_FLUSH);
		xdrs->x_op = XDR_FREE;
		xdr_bytes(xdrs, (char **) &zbuf, &size,
				objp->compression.cblksize);
		xdrs->x_op = XDR_DECODE;
		if (status != Z_OK) {
			break;
		}
	} while (size == objp->compression.cblksize);

	do {
		status = inflate(&z, Z_FINISH);
	} while (status == Z_OK);

	/*
	 * Now, there should be a decompressed data stream in obuf.
	 * Run the standard xdr routines using a memory xdr now.
	 */

	inflateEnd(&z);
	free(zbuf);
	xdr_destroy(&mem_xdrs);

	xdrmem_create(&mem_xdrs, obuf, objp->compression.cfullsize, XDR_DECODE);
	status = xdr_ecs_Result_Work(&mem_xdrs, objp);
	xdr_destroy(&mem_xdrs);

	return status;
}


bool_t
xdr_ecs_Result_Encode(XDR *xdrs, ecs_Result *objp)
{
	XDR mem_xdrs;
	int attempts;
	void *zbuf;
	z_stream z;
	int count;
	int status;
	
	if (objp->compression.ctype == ECS_COMPRESS_NONE) {
		if (!xdr_u_int(xdrs, &objp->compression.cfullsize)) {
			return FALSE;
		}
		return xdr_ecs_Result_Work(xdrs, objp);
	}
	if (objp->compression.ctype != ECS_COMPRESS_ZLIB) {
		return (FALSE);
	}
	if (objp->compression.cblksize == 0) {
		if (!xdr_u_int(xdrs, &objp->compression.cfullsize)) {
			return FALSE;
		}
		return xdr_ecs_Result_Work(xdrs, objp);
	}
	attempts = 0;
	while (attempts < 7) {
		if (enough == 0) {
			obufsize = obufsize * 2 + 200000;
			if (obuf) free(obuf);
			obuf = malloc(obufsize);
			if (obuf == NULL) {
				obufsize = 0;
				return (FALSE);
			}
		}
		if (attempts > 0) {
			xdr_destroy(&mem_xdrs);
		}
		xdrmem_create(&mem_xdrs, obuf, obufsize, XDR_ENCODE);
	
		enough = xdr_ecs_Result_Work(&mem_xdrs, objp);
		if (enough) {
			break;
		}
		attempts++;
	}
	if (! enough) {
		xdr_destroy(&mem_xdrs);
		free(obuf);
		obuf = NULL;
		return (FALSE);
	}
	zbuf = malloc(objp->compression.cblksize);
	if (zbuf == NULL) {
		xdr_destroy(&mem_xdrs);
		return (FALSE);
	}
	
	objp->compression.cfullsize = XDR_GETPOS(&mem_xdrs);
	if (!xdr_u_int(xdrs, &objp->compression.cfullsize)) {
		xdr_destroy(&mem_xdrs);
		return FALSE;
	}

	/*
	 * Basic algorithm is this: call the compress routine until
	 * the output buffer has filled up with compression.cblksize
	 * bytes.  Then write these bytes out with XDR_PUT_BYTES().
	 * Continue doing this until there are no more bytes left.
	 * The final block must be less than compression.cblksize bytes.
	 * If it is exactly equal to compression.cblksize bytes, then
	 * a final 0 length block needs to be sent with XDR_PUT_BYTES().
	 */
	z.zalloc = NULL;
	z.zfree = NULL;
	z.opaque = NULL;
	if (deflateInit(&z, objp->compression.clevel) != Z_OK) {
		xdr_destroy(&mem_xdrs);
		free(zbuf);
		return (FALSE);
	}
	
	z.next_in = obuf;
	z.avail_in = objp->compression.cfullsize;
	while (1) {
		z.next_out = zbuf;
		z.avail_out = objp->compression.cblksize;
		status = deflate(&z, Z_NO_FLUSH);
		count = objp->compression.cblksize - z.avail_out;
		if ((status != Z_OK) || (count < objp->compression.cblksize)) {
			break;
		}
		xdr_bytes(xdrs, (char **) &zbuf, &count,
			  objp->compression.cblksize);
	}
	while (1) {
		status = deflate(&z, Z_FINISH);
		count = objp->compression.cblksize - z.avail_out;
		if ((status != Z_OK) || (count < objp->compression.cblksize)) {
			xdr_bytes(xdrs, (char **) &zbuf, &count,
				  objp->compression.cblksize);
			if (count == objp->compression.cblksize) {
				count = 0;
				xdr_bytes(xdrs, (char **) &zbuf, &count,
					  objp->compression.cblksize);
			}
			break;
		}
		xdr_bytes(xdrs, (char **) &zbuf, &count,
			  objp->compression.cblksize);
		z.next_out = zbuf;
		z.avail_out = objp->compression.cblksize;
	}
	deflateEnd(&z);
	free(zbuf);
	xdr_destroy(&mem_xdrs);

	return (TRUE);
}

bool_t
xdr_ecs_Result(XDR *xdrs, ecs_Result *objp)
{
	 if (!xdr_u_int(xdrs, &objp->compression.cachesize)) {
		 return FALSE;
	 }
	 if (!xdr_u_int(xdrs, &objp->compression.ctype)) {
		 return FALSE;
	 }
	 if (!xdr_u_int(xdrs, &objp->compression.cversion)) {
		 return FALSE;
	 }
	 if (!xdr_u_int(xdrs, &objp->compression.cblksize)) {
		 return FALSE;
	 }
	 if (xdrs->x_op == XDR_ENCODE) {
		 return xdr_ecs_Result_Encode(xdrs, objp);
	 } else if (xdrs->x_op == XDR_DECODE) {
		 return xdr_ecs_Result_Decode(xdrs, objp);
	 } else {
		 return xdr_ecs_Result_Free(xdrs, objp);
	 }
}

