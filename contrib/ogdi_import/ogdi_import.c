#include "ecs.h"
#include "shapefil.h"
#include <assert.h>


static int	ClientID = -1;

/************************************************************************/
/*                             CheckError()                             */
/************************************************************************/

static int CheckError( ecs_Result * result )

{
    if( ECSERROR( result ) ) {
        printf( "ECSERROR: %s\n", result->message );
        return TRUE;
    }
    else
        return FALSE;
}

/************************************************************************/
/*                             AccessURL()                              */
/************************************************************************/

static int AccessURL( char * url, ecs_Region * region )

{
    ecs_Result *result;

/* -------------------------------------------------------------------- */
/*      Close old client if there is one active.                        */
/* -------------------------------------------------------------------- */
    if( ClientID != -1 ) {
        result = cln_DestroyClient(ClientID);
        if( CheckError( result ) )
            return( FALSE );
        ClientID = -1;
    }
            
/* -------------------------------------------------------------------- */
/*      Open new URL/client.                                            */
/* -------------------------------------------------------------------- */
    result = cln_CreateClient(&ClientID,url);
    if( CheckError( result ) )
        return( FALSE );

/* -------------------------------------------------------------------- */
/*      print the global bounds.                                        */
/* -------------------------------------------------------------------- */
    result = cln_GetGlobalBound(ClientID);
    if( CheckError( result ) )
        return( FALSE );

    printf( "Global Bounds\n" );
    printf( "north = %f\n", ECSREGION(result).north );
    printf( "south = %f\n", ECSREGION(result).south );
    printf( "east = %f\n", ECSREGION(result).east );
    printf( "west = %f\n", ECSREGION(result).west );
    printf( "ns_res = %f\n", ECSREGION(result).ns_res );
    printf( "ew_res = %f\n", ECSREGION(result).ew_res );

    *region = ECSREGION(result);

/* -------------------------------------------------------------------- */
/*      Print the Dictionary (update).                                  */
/* -------------------------------------------------------------------- */
    result = cln_UpdateDictionary( ClientID, "" );
    if( CheckError( result ) )
        return( FALSE );

    printf( "UpdateDictionary = \n%s\n", ECSTEXT(result) );
  
/* -------------------------------------------------------------------- */
/*      Print the projection.                                           */
/* -------------------------------------------------------------------- */
    result = cln_GetServerProjection( ClientID );
    if( CheckError( result ) )
        return( FALSE );

    printf( "Projection = `%s'\n", ECSTEXT(result) );
  
    return( TRUE );
}

/************************************************************************/
/*                            ImportVectors()                           */
/************************************************************************/

static void ImportVectors( ecs_Region *region, const char * layer,
                           const char * out_file,
                           ecs_Family featureType )

{
    ecs_Result *result;
    ecs_LayerSelection selectionLayer;
    SHPHandle   hShape;
    DBFHandle   hDBF;
    char	filename[1024];
    ecs_ObjAttributeFormat *oaf;
    int		i, field_count;

/* -------------------------------------------------------------------- */
/*      Select a region ... this should be overridable from the         */
/*      command line.                                                   */
/* -------------------------------------------------------------------- */
    result = cln_SelectRegion(ClientID,region);
    if( CheckError( result ) )
        return;

/* -------------------------------------------------------------------- */
/*      Define the layer to select.  For now we only support lines.     */
/* -------------------------------------------------------------------- */
    assert( featureType == Line );
    
    selectionLayer.Select = (char *) layer;
    selectionLayer.F = featureType;
    
    result = cln_SelectLayer(ClientID,&selectionLayer);
    
    if( CheckError(result) )
        return;

/* -------------------------------------------------------------------- */
/*      Create the shapefile to write the lines to.                     */
/* -------------------------------------------------------------------- */
    sprintf( filename, "%s.shp", out_file );
    hShape = SHPCreate( filename, SHPT_ARC );

    sprintf( filename, "%s.dbf", out_file );
    hDBF = DBFCreate( filename );

/* -------------------------------------------------------------------- */
/*      Setup DBF schema.                                               */
/* -------------------------------------------------------------------- */
    result = cln_GetAttributesFormat( ClientID );
    if( CheckError( result ) )
        return;
        
    oaf = &(ECSRESULT(result).oaf);
    field_count = oaf->oa.oa_len;
    for( i = 0; i < field_count; i++ )
    {
        switch( oaf->oa.oa_val[i].type )
        {
          case Char:
          case Varchar:
          case Longvarchar:
            if( oaf->oa.oa_val[i].lenght > 0 )
                DBFAddField( hDBF, oaf->oa.oa_val[i].name,
                             FTString, oaf->oa.oa_val[i].lenght, 0 );
            else
                DBFAddField( hDBF, oaf->oa.oa_val[i].name,
                             FTString, 64, 0 );
            break;

          case Decimal:
          case Smallint:
          case Integer:
            if( oaf->oa.oa_val[i].lenght > 0 )
                DBFAddField( hDBF, oaf->oa.oa_val[i].name,
                             FTDouble, oaf->oa.oa_val[i].lenght, 0 );
            else 
                DBFAddField( hDBF, oaf->oa.oa_val[i].name,
                             FTString, 11, 0 );
            break;

          case Numeric:
          case Real:
          case Float:
          case Double:
            if( oaf->oa.oa_val[i].lenght > 0 )
                DBFAddField( hDBF, oaf->oa.oa_val[i].name,
                             FTDouble,
                             oaf->oa.oa_val[i].lenght,
                             oaf->oa.oa_val[i].precision );
            else
                DBFAddField( hDBF, oaf->oa.oa_val[i].name,
                             FTString, 18, 7 );
            break;

        }
    }
    

/* -------------------------------------------------------------------- */
/*      Process objects.                                                */
/* -------------------------------------------------------------------- */
    result = cln_GetNextObject(ClientID);
    while (ECSSUCCESS(result))
    {
        char	*pszList;
        int	iField, iShape;
        
/* -------------------------------------------------------------------- */
/*      Write shapefile geometry                                        */
/* -------------------------------------------------------------------- */
        
        if( featureType == Line ) {
            ecs_Line	*line = &(ECSGEOM(result).line);
            double	*x, *y;
            SHPObject   *psShape;

            x = (double *) malloc(sizeof(double) * line->c.c_len);
            y = (double *) malloc(sizeof(double) * line->c.c_len);

            for( i=0; i < (int) line->c.c_len; i++ ) {
                x[i] = line->c.c_val[i].x;
                y[i] = line->c.c_val[i].y;
            }

            psShape = SHPCreateSimpleObject( SHPT_ARC, line->c.c_len,
                                             x, y, NULL );

            free( x );
            free( y );
            
            iShape = SHPWriteObject( hShape, -1, psShape );

            SHPDestroyObject( psShape );
        }

/* -------------------------------------------------------------------- */
/*      Write attributes to dbf file.                                   */
/* -------------------------------------------------------------------- */
        pszList = ECSOBJECTATTR(result);
        
        for( iField = 0; iField < field_count; iField++ )
        {
            char	*pszFieldStart;
            int		nNameLen;
            char	chSavedChar;
            
            /* parse out the next attribute value */
            if( !ecs_FindElement( pszList, &pszFieldStart, &pszList,
                                  &nNameLen, NULL ) )
            {
                nNameLen = 0;
                pszFieldStart = pszList;
            }

            /* Skip any trailing white space (for string constants). */

            if( nNameLen > 0 && pszFieldStart[nNameLen-1] == ' ' )
                nNameLen--;

	    /* zero terminate the single field value, but save the           */
            /* character we overwrote, so we can restore it when done.       */
            
            chSavedChar = pszFieldStart[nNameLen];
            pszFieldStart[nNameLen] = '\0';

            switch( DBFGetFieldInfo( hDBF, iField, NULL, NULL, NULL ) )
            {
              case FTString:
                DBFWriteStringAttribute( hDBF, iShape, iField, pszFieldStart );
                break;

              case FTInteger:
                DBFWriteIntegerAttribute( hDBF, iShape, iField,
                                          atoi(pszFieldStart) );
                break;

              case FTDouble:
                DBFWriteIntegerAttribute( hDBF, iShape, iField,
                                          atof(pszFieldStart) );
                break;

              default:
		break;
            }

            pszFieldStart[nNameLen] = chSavedChar;
        }
        
/* -------------------------------------------------------------------- */
/*      Read next object                                                */
/* -------------------------------------------------------------------- */
        result = cln_GetNextObject(ClientID);
    }

    SHPClose( hShape );
    DBFClose( hDBF );
}

/************************************************************************/
/*                            ImportMatrix()                            */
/************************************************************************/

static void ImportMatrix( ecs_Region *region, const char * layer,
                          const char * out_file )

{
    FILE	*fp_aux, *fp_raw;
    char	filename[256];
    int		i, xsize, ysize;
    ecs_Result *result;
    ecs_LayerSelection selectionLayer;

/* -------------------------------------------------------------------- */
/*      Select a region ... this should be overridable from the         */
/*      command line.                                                   */
/* -------------------------------------------------------------------- */
    result = cln_SelectRegion(ClientID,region);
    if( CheckError( result ) )
        return;

    xsize = (region->east - region->west) / region->ew_res;
    ysize = (region->north - region->south) / region->ns_res;

/* -------------------------------------------------------------------- */
/*      Define the layer to select.                                     */
/* -------------------------------------------------------------------- */
    selectionLayer.Select = (char *) layer;
    selectionLayer.F = Matrix;
    
    result = cln_SelectLayer(ClientID,&selectionLayer);
    
    if( CheckError( result ) )
        return;

/* -------------------------------------------------------------------- */
/*	Write .aux file header.						*/
/* -------------------------------------------------------------------- */ 
    sprintf( filename, "%s.aux", out_file );
    fp_aux = fopen( filename, "wt" );

    sprintf( filename, "%s.raw", out_file );

    fprintf( fp_aux, "AuxilaryTarget: %s\n", filename );

/* -------------------------------------------------------------------- */
/*      Get the raster information.                                     */
/* -------------------------------------------------------------------- */
    result = cln_GetRasterInfo( ClientID );
    if( CheckError( result ) )
        return;

    fprintf( fp_aux, "RawDefinition: %d %d 1\n", xsize, ysize );

    fprintf( fp_aux, "ChanDefinition-1: 8U 0 1 %d Swapped\n", xsize );

    for( i = 0; i < (int) ECSRASTERINFO(result).cat.cat_len; i++ ) {
        fprintf( fp_aux, "METADATA_IMG_1_Class_%d_Name: %s\n",
                 (int) ECSRASTERINFO(result).mincat + i, 
                 ECSRASTERINFO(result).cat.cat_val[i].label );
        fprintf( fp_aux, "METADATA_IMG_1_Class_%d_Color: (RGB:%d %d %d)\n",
                 (int) ECSRASTERINFO(result).mincat + i, 
                 ECSRASTERINFO(result).cat.cat_val[i].r,
                 ECSRASTERINFO(result).cat.cat_val[i].g,
                 ECSRASTERINFO(result).cat.cat_val[i].b );
    }

/* -------------------------------------------------------------------- */
/*      Attach the projection information.                              */
/* -------------------------------------------------------------------- */
    fprintf( fp_aux, "MapUnits: METER\n" );
    fprintf( fp_aux, "UpLeftX: %24.15E\n", region->west );
    fprintf( fp_aux, "UpLeftY: %25.15E\n", region->north );
    fprintf( fp_aux, "LoRightX: %25.15E\n", region->east );
    fprintf( fp_aux, "LoRightY: %25.15E\n", region->south );

    fclose( fp_aux );

/* -------------------------------------------------------------------- */
/*      open the raw data file.                                         */
/* -------------------------------------------------------------------- */
    fp_raw = fopen( filename, "wb" );
    
    result = cln_GetNextObject(ClientID);
    while (ECSSUCCESS(result)) {
        
        static int report_flag = 0;
        
        if( (int) ECSOBJECT(result).geom.ecs_Geometry_u.matrix.x.x_len
            != xsize && report_flag == 0 ) {

            report_flag = 1;
            printf( "Got %d pixels, instead of expected %d pixels!\n",
                    (int) ECSOBJECT(result).geom.ecs_Geometry_u.matrix.x.x_len,
                    xsize );
        }
            
        for( i = 0; i < xsize; i++ )
        {
            fprintf( fp_raw, "%c", ECSRASTER(result)[i] );
        }

        result = cln_GetNextObject(ClientID);
    }

    fclose( fp_raw );
}

/************************************************************************/
/*                            ImportImage()                             */
/************************************************************************/

static void ImportImage( ecs_Region *region, const char * layer,
                         const char * out_file )

{
    FILE	*fp_aux, *fp_raw;
    char	filename[256];
    int		i, xsize, ysize, width_code;
    ecs_Result *result;
    ecs_LayerSelection selectionLayer;

/* -------------------------------------------------------------------- */
/*      Select a region ... this should be overridable from the         */
/*      command line.                                                   */
/* -------------------------------------------------------------------- */
    result = cln_SelectRegion(ClientID,region);
    if( CheckError( result ) )
        return;

    xsize = (region->east - region->west) / region->ew_res;
    ysize = (region->north - region->south) / region->ns_res;

/* -------------------------------------------------------------------- */
/*      Define the layer to select.                                     */
/* -------------------------------------------------------------------- */
    selectionLayer.Select = (char *) layer;
    selectionLayer.F = Image;
    
    result = cln_SelectLayer(ClientID,&selectionLayer);
    
    if( CheckError( result ) )
        return;

/* -------------------------------------------------------------------- */
/*	Write .aux file header.						*/
/* -------------------------------------------------------------------- */ 
    sprintf( filename, "%s.aux", out_file );
    fp_aux = fopen( filename, "wt" );

    sprintf( filename, "%s.raw", out_file );

    fprintf( fp_aux, "AuxilaryTarget: %s\n", filename );

/* -------------------------------------------------------------------- */
/*      Get the raster information.                                     */
/* -------------------------------------------------------------------- */
    result = cln_GetRasterInfo( ClientID );
    if( CheckError( result ) )
        return;

    width_code = ECSRASTERINFO(result).width;
    if( width_code == 1 /* RGBA */ ) {
        fprintf( fp_aux, "RawDefinition: %d %d 4\n", xsize, ysize );
        fprintf( fp_aux, "ChanDefinition-1: 8U 0 4 %d\n", xsize*4 );
        fprintf( fp_aux, "ChanDefinition-2: 8U 1 4 %d\n", xsize*4 );
        fprintf( fp_aux, "ChanDefinition-3: 8U 2 4 %d\n", xsize*4 );
        fprintf( fp_aux, "ChanDefinition-4: 8U 3 4 %d\n", xsize*4 );
    } else if( width_code == 5 /* int32 */ ) {
        fprintf( fp_aux, "RawDefinition: %d %d 1\n", xsize, ysize );
        fprintf( fp_aux, "ChanDefinition-1: 16U 0 2 %d Swapped\n", xsize*2 );
    }
    else {
        printf( "Unsupported Image raster type %d.\n",
                ECSRASTERINFO(result).width );
        assert( FALSE );
        return;
    }

/* -------------------------------------------------------------------- */
/*      Attach the projection information.                              */
/* -------------------------------------------------------------------- */
    fprintf( fp_aux, "MapUnits: METER\n" );
    fprintf( fp_aux, "UpLeftX: %24.15E\n", region->west );
    fprintf( fp_aux, "UpLeftY: %25.15E\n", region->north );
    fprintf( fp_aux, "LoRightX: %25.15E\n", region->east );
    fprintf( fp_aux, "LoRightY: %25.15E\n", region->south );

    fclose( fp_aux );

/* -------------------------------------------------------------------- */
/*      open the raw data file.                                         */
/* -------------------------------------------------------------------- */
    fp_raw = fopen( filename, "wb" );
    
    result = cln_GetNextObject(ClientID);
    while (ECSSUCCESS(result)) {
        
        static int report_flag = 0;
        
        if( (int) ECSOBJECT(result).geom.ecs_Geometry_u.matrix.x.x_len
            != xsize && report_flag == 0 ) {

            report_flag = 1;
            printf( "Got %d pixels, instead of expected %d pixels!\n",
                    (int) ECSOBJECT(result).geom.ecs_Geometry_u.matrix.x.x_len,
                    xsize );
        }
            
        for( i = 0; i < xsize; i++ ) {
            if( width_code == 5 ) {
                fprintf( fp_raw, "%c%c",
                         (ECSRASTER(result)[i] & 0xff),
                         (ECSRASTER(result)[i] >> 8) & 0xff );
            } else {
#ifdef notdef
                unsigned char red, green, blue, trans;
                
                ogdi_GetRGBFromPixel( ECSRASTER(result)[i],
                                      &trans, &red, &green, &blue );
                fprintf( fp_raw, "%c%c%c%c",
                         red, green, blue, trans );
#endif                
            }
        }

        result = cln_GetNextObject(ClientID);
    }

    fclose( fp_raw );
}

/************************************************************************/
/*                                main()                                */
/************************************************************************/

int main( int argc, char ** argv )
{
    ecs_Family	featureType = Point;
    char *layer = "";
    static ecs_Region	reg, *region;
    ecs_Result *result;
    int		i;
    char	*out_file = "ogdi_out";

    if( argc == 1 )
    {
        printf("Usage: ogdi_import -u url -f family\n");
        printf("          [-r north south east west] [-res ns_res ew_res]\n" );
        printf("          [-o filename]\n" );
        printf("          -l layername [more_opts -l layername]\n" );
        exit( 1 );
    }

/* -------------------------------------------------------------------- */
/*      Handle commandline arguments.                                   */
/* -------------------------------------------------------------------- */
    for( i = 1; i < argc; i++ ) {

        if( i == argc - 1 ) {
            /* skip ... the rest require arguments.  */
        }
        else if( strcmp(argv[i],"-u") == 0 ) {
            AccessURL( argv[++i], &reg );
            if( region == NULL )
                region = &reg;
        }
        else if( strcmp(argv[i],"-o") == 0 ) {
            out_file = argv[++i];
        }
        else if( strcmp(argv[i], "-l") == 0 ) {
            layer = argv[++i];

            if( featureType == Matrix ) {
                ImportMatrix( region, layer, out_file );
            }
            else if( featureType == Image ) {
                ImportImage( region, layer, out_file );
            }
            else {
                ImportVectors( region, layer, out_file, featureType );
            }
        }
        else if( strcmp(argv[i], "-f") == 0 ) {
            if( strcmp(argv[i+1],"Point") == 0 )
                featureType = Point;
            else if( strcmp(argv[i+1],"Line") == 0 )
                featureType = Line;
            else if( strcmp(argv[i+1],"Area") == 0 )
                featureType = Area;
            else if( strcmp(argv[i+1],"Text") == 0 )
                featureType = Text;
            else if( strcmp(argv[i+1],"Matrix") == 0 )
                featureType = Matrix;
            else if( strcmp(argv[i+1],"Image") == 0 )
                featureType = Image;
            i++;
        }
        else if( strcmp(argv[i], "-r") == 0 && i < argc - 4 ) {

            region = &reg;
            reg.north = atof(argv[++i]);
            reg.south = atof(argv[++i]);
            reg.east = atof(argv[++i]);
            reg.west = atof(argv[++i]);
        }
        else if( strcmp(argv[i], "-res") == 0 && i < argc - 2 ) {

            reg.ns_res = atof(argv[++i]);
            reg.ew_res = atof(argv[++i]);
        }
    }

/* -------------------------------------------------------------------- */
/*      Close old client if there is one active.                        */
/* -------------------------------------------------------------------- */
    if( ClientID != -1 ) {
        result = cln_DestroyClient(ClientID);
        if( CheckError( result ) )
            return( FALSE );
        ClientID = -1;
    }

    return 0;
}

