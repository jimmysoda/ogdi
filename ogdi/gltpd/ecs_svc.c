/*
 * Please do not edit this file.
 * It was generated using rpcgen.
 */

#include "ecs.h"
#ifndef _WINDOWS
#include <stdio.h>
#include <stdlib.h>/* getenv, exit */
#include <rpc/pmap_clnt.h> /* for pmap_unset */
#include <string.h> /* strcmp */ 
#include <signal.h>
#include <sys/ioctl.h> /* ioctl, TIOCNOTTY */
#include <sys/types.h> /* open */
#include <sys/stat.h> /* open */
#include <fcntl.h> /* open */
#include <unistd.h> /* getdtablesize */
#include <memory.h>
#include <sys/socket.h>
#include <netinet/in.h>
#endif /* ifndef _WINDOWS */

extern unsigned long newprogramno;
#include <syslog.h>

#ifndef SIG_PF
#define	SIG_PF void(*)(int)
#endif

#define _RPCSVC_CLOSEDOWN 120
#define main dummy
#define RPC_SVC_FG

extern ecs_Server *svr_handle;

static void closedown(int sig);

void start_closedown_check(void)
{
#ifndef _WINDOWS
	(void) signal(SIGALRM, (SIG_PF) closedown);
	(void) alarm(_RPCSVC_CLOSEDOWN);
#endif /* ifndef _WINDOWS */
}

 int _rpcpmstart;		/* Started by a port monitor ? */
 int _rpcfdtype;		/* Whether Stream or Datagram ? */
 int _rpcsvcdirty;	/* Still serving ? */

static
void _msgout(char* msg)
{
#ifdef RPC_SVC_FG
	if (_rpcpmstart)
		(void) fprintf(stderr, "%s\n",  msg);
	else
		(void) fprintf(stderr, "%s\n", msg);
#else
	(void) fprintf(stderr, "%s\n",  msg);
#endif
}

static void
closedown(int sig)
{
	(void) signal(sig, (SIG_PF) closedown);
	if (_rpcsvcdirty == 0) {
		extern fd_set svc_fdset;
		static int size;
		int i, openfd;

		if (_rpcfdtype == SOCK_DGRAM)
			exit(0);
		if (size == 0) {
			size = getdtablesize();
		}
		for (i = 0, openfd = 0; i < size && openfd < 2; i++)
			if (FD_ISSET(i, &svc_fdset))
				openfd++;
		if (openfd <= 1)
			exit(0);
	}
	(void) alarm(_RPCSVC_CLOSEDOWN);
}

void
ecsprog_1(struct svc_req *rqstp, register SVCXPRT *transp)
{
	union {
		char *createserver_1_arg;
		ecs_LayerSelection selectlayer_1_arg;
		ecs_LayerSelection releaselayer_1_arg;
		ecs_Region selectregion_1_arg;
		char *getobject_1_arg;
		ecs_Coordinate getobjectidfromcoord_1_arg;
		char *updatedictionary_1_arg;
		u_int setserverlanguage_1_arg;
		char *setserverprojection_1_arg;
		ecs_RasterConversion setrasterconversion_1_arg;
		ecs_ProxyCreateServer createproxyserver_1_arg;
		ecs_Compression setcompression_1_arg;
	} argument;
	char *result;
	xdrproc_t xdr_argument, xdr_result;
	char *(*local)(char *, struct svc_req *);

	_rpcsvcdirty = 1;
	switch (rqstp->rq_proc) {
	case NULLPROC:
		(void) svc_sendreply(transp, (xdrproc_t) xdr_void, (char *)NULL);
		_rpcsvcdirty = 0;
		return;

	case CREATESERVER:
		xdr_argument = (xdrproc_t) xdr_wrapstring;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) createserver_1_svc;
		break;

	case DESTROYSERVER:
		xdr_argument = (xdrproc_t) xdr_void;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) destroyserver_1_svc;
		break;

	case SELECTLAYER:
		xdr_argument = (xdrproc_t) xdr_ecs_LayerSelection;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) selectlayer_1_svc;
		break;

	case RELEASELAYER:
		xdr_argument = (xdrproc_t) xdr_ecs_LayerSelection;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) releaselayer_1_svc;
		break;

	case SELECTREGION:
		xdr_argument = (xdrproc_t) xdr_ecs_Region;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) selectregion_1_svc;
		break;

	case GETDICTIONNARY:
		xdr_argument = (xdrproc_t) xdr_void;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getdictionnary_1_svc;
		break;

	case GETATTRIBUTEFORMAT:
		xdr_argument = (xdrproc_t) xdr_void;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getattributeformat_1_svc;
		break;

	case GETNEXTOBJECT:
		xdr_argument = (xdrproc_t) xdr_void;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getnextobject_1_svc;
		break;

	case GETRASTERINFO:
		xdr_argument = (xdrproc_t) xdr_void;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getrasterinfo_1_svc;
		break;

	case GETOBJECT:
		xdr_argument = (xdrproc_t) xdr_wrapstring;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getobject_1_svc;
		break;

	case GETOBJECTIDFROMCOORD:
		xdr_argument = (xdrproc_t) xdr_ecs_Coordinate;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getobjectidfromcoord_1_svc;
		break;

	case UPDATEDICTIONARY:
		xdr_argument = (xdrproc_t) xdr_wrapstring;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) updatedictionary_1_svc;
		break;

	case GETSERVERPROJECTION:
		xdr_argument = (xdrproc_t) xdr_void;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getserverprojection_1_svc;
		break;

	case GETGLOBALBOUND:
		xdr_argument = (xdrproc_t) xdr_void;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) getglobalbound_1_svc;
		break;

	case SETSERVERLANGUAGE:
		xdr_argument = (xdrproc_t) xdr_u_int;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) setserverlanguage_1_svc;
		break;

	case SETSERVERPROJECTION:
		xdr_argument = (xdrproc_t) xdr_wrapstring;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) setserverprojection_1_svc;
		break;

	case SETRASTERCONVERSION:
		xdr_argument = (xdrproc_t) xdr_ecs_RasterConversion;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) setrasterconversion_1_svc;
		break;

	case CREATEPROXYSERVER:
		xdr_argument = (xdrproc_t) xdr_ecs_ProxyCreateServer;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) createproxyserver_1_svc;
		break;

	case SETCOMPRESSION:
		xdr_argument = (xdrproc_t) xdr_ecs_Compression;
		xdr_result = (xdrproc_t) xdr_ecs_Result;
		local = (char *(*)(char *, struct svc_req *)) setcompression_1_svc;
		break;

	default:
		svcerr_noproc(transp);
		_rpcsvcdirty = 0;
		return;
	}
	(void) memset((char *)&argument, 0, sizeof (argument));
	if (!svc_getargs(transp, xdr_argument, (caddr_t) &argument)) {
		svcerr_decode(transp);
		_rpcsvcdirty = 0;
		return;
	}
	result = (*local)((char *)&argument, rqstp);
	if (result) {
	  ecs_Result *tmp = (ecs_Result *) result;
	  /* The CreateServer call must be always uncompressed */
	  if (svr_handle && rqstp->rq_proc != CREATESERVER) {
	    tmp->compression = svr_handle->compression;
	  } else {
	    /* This will occur after a DESTROYSERVER call */
	    tmp->compression.ctype = ECS_COMPRESS_NONE;
	    tmp->compression.cversion = 0;
	    tmp->compression.clevel = 0;
	    tmp->compression.cblksize = 0;
	    tmp->compression.cfullsize = 0;
	    tmp->compression.cachesize = 0;
	  }		
	}
	if (result != NULL && !svc_sendreply(transp, xdr_result, result)) {
		svcerr_systemerr(transp);
	}
	if (!svc_freeargs(transp, xdr_argument, (caddr_t) &argument)) {
		_msgout("unable to free arguments");
		exit(1);
	}

	if (rqstp->rq_proc == DESTROYSERVER) {
#ifdef _WINDOWS
	  rpc_nt_exit();
#endif
	  pmap_unset(newprogramno, ECSVERS);
	  exit(0);
	}

	_rpcsvcdirty = 0;
	return;
}

