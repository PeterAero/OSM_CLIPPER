#include <QCoreApplication>
#include <sqlite3.h>
#include <iostream>
#include <string>
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogrsf_frmts.h"
#include <ogr_geometry.h>



void computeAngle(double& angle, OGRPoint& Point1, OGRPoint& Point2){

    angle = std::atan2(Point1.getY() - Point2.getY(), Point1.getX() - Point2.getX());

}

void loadCompleteOSMToMem(OGRPolygon& MyRing){

    GDALDataset * poDS;
    poDS = (GDALDataset*) GDALOpenEx( "/home/peter/Downloads/hessen/gis.osm_roads_free_1.shp", GDAL_OF_VECTOR, NULL, NULL, NULL );

    if( poDS == NULL )
    {
        printf( "Open failed.\n" );
        exit( 1 );
    }

    std::cout << poDS->GetLayerCount() << std::endl;

    OGRGeometry * NewGeometry;
    NewGeometry = &MyRing;

    OGRLayer * poLayer;
    poLayer = poDS->GetLayer(0);
    NewGeometry->transformTo(poLayer->GetSpatialRef());

    poLayer->SetSpatialFilter(NewGeometry);

    OGRFeature *poFeature;
    poLayer->ResetReading();

    // * * * * * * * * * * * * * * * * print out db content * * * * * * * * * * * * * * * *
    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {

        OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
        int iField;
        for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ )
        {
            OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
            if( poFieldDefn->GetType() == OFTInteger )
                printf( "%d,", poFeature->GetFieldAsInteger( iField ) );
            else if( poFieldDefn->GetType() == OFTInteger64 )
                printf( CPL_FRMT_GIB ",", poFeature->GetFieldAsInteger64( iField ) );
            else if( poFieldDefn->GetType() == OFTReal )
                printf( "%.3f,", poFeature->GetFieldAsDouble(iField) );
            else if( poFieldDefn->GetType() == OFTString )
                printf( "%s,", poFeature->GetFieldAsString(iField) );
            else
                printf( "%s,", poFeature->GetFieldAsString(iField) );
        }
        OGRGeometry *poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL
                && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )
        {
            OGRPoint *poPoint = (OGRPoint *) poGeometry;
            printf( "%.3f,%3.f\n", poPoint->getX(), poPoint->getY() );
        }
        else
        {
            printf( "no point geometry\n" );
        }
        OGRFeature::DestroyFeature( poFeature );
    }

    GDALClose( poDS );

}


void createBoundingBox(const std::string& AbsPath, OGRLinearRing& MyRing){

    /* function to create the bounding box of a given image
     * as OGRLinearRing object in geographic coordinates
     *
     * */

    // * * * * * get corner coordinates of image in image coordinate system * * * * * * *
    GDALDataset  * srcDataset;
    srcDataset = (GDALDataset *) GDALOpen( AbsPath.c_str(), GA_ReadOnly );
    double geoParam[6];
    double x[4];
    double y[4];

    srcDataset->GetGeoTransform(geoParam);
    double Cols = srcDataset->GetRasterXSize();
    double Rows = srcDataset->GetRasterYSize();

    x[0] = geoParam[0];
    x[1] = geoParam[0] + Cols * geoParam[1];
    x[2] = geoParam[0] + Cols * geoParam[1];
    x[3] = geoParam[0];
    y[0] = geoParam[3];
    y[1] = geoParam[3];
    y[2] = geoParam[3] + Rows * geoParam[5];
    y[3] = geoParam[3] + Rows * geoParam[5];

    // * * * * * * * * * * * * * * create SRS from input image * * * * * * * * * * * * * *
    OGRSpatialReference oSourceSRS;
    const char * poWKT = srcDataset->GetProjectionRef();
    char * poWKT_tmp = (char *) poWKT;
    oSourceSRS.importFromWkt(&poWKT_tmp);

    // * * * * * * * create transformation object from input image to polygon SRS * * * * * * *
    OGRCoordinateTransformation *poCT;
    poCT = OGRCreateCoordinateTransformation( &oSourceSRS, MyRing.getSpatialReference() );
    poCT->Transform(4, x, y);

    OGRPoint MyPoint1(x[0], y[0]);
    OGRPoint MyPoint2(x[1], y[1]);
    OGRPoint MyPoint3(x[2], y[2]);
    OGRPoint MyPoint4(x[3], y[3]);

    MyRing.addPoint(&MyPoint1);
    MyRing.addPoint(&MyPoint2);
    MyRing.addPoint(&MyPoint3);
    MyRing.addPoint(&MyPoint4);
    MyRing.addPoint(&MyPoint1);

    GDALClose(srcDataset);

}

void createTmpLayer(OGRPolygon& MyRing, const std::string& tmpPath_Rectangle){

    /* function which creates a shapefile consisting of a single rectangle,
     * usually the foodprint of an aerial/satellite imagery
     *
     * */

    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ){
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }

    GDALDataset *poDS;
    poDS = poDriver->Create( tmpPath_Rectangle.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ){
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer( "point_out", MyRing.getSpatialReference(), wkbPolygon, NULL );

    if( poLayer == NULL ){
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    OGRFieldDefn oField( "Name", OFTString );
    oField.SetWidth(32);
    if( poLayer->CreateField( &oField ) != OGRERR_NONE ){
        printf( "Creating Name field failed.\n" );
        exit( 1 );
    }

    OGRFeature *poFeature;
    poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
    poFeature->SetField( "Name", "whatever" );

    poFeature->SetGeometry( &MyRing );
    if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ){
        printf( "Failed to create feature in shapefile.\n" );
        exit( 1 );
    }
    OGRFeature::DestroyFeature( poFeature );

    GDALClose( poDS );
}

void clipDB(const std::string& AbsPathDB, const std::string& tmpPath_Rectangle,
            const std::string& tmpPath_clip){

    /* function to clip the data of a given vector file (e.g. OSM DB) with respect to
     * a second vector file (e.g. footprint of an image)
     *
     *
     * */

    // * * * * * * * * * * * * * * * * * * load OSM DB * * * * * * * * * * * * * * * * * *
    GDALDataset * poDS1;
    poDS1 = (GDALDataset*) GDALOpenEx( AbsPathDB.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );

    if( poDS1 == NULL ){
        printf( "Open of DB failed.\n" );
        exit( 1 );
    }
    OGRLayer * poDBLayer;
    poDBLayer = poDS1->GetLayer(0);

    // * * * * * * * * * * * * * * * * * * load Mask Layer * * * * * * * * * * * * * * * * * *
    GDALDataset * poDS2;
    poDS2 = (GDALDataset*) GDALOpenEx( tmpPath_Rectangle.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );
    if( poDS2 == NULL ){
        printf( "Open of Rectangle file failed.\n" );
        exit( 1 );
    }
    OGRLayer * poMaskLayer;
    poMaskLayer = poDS2->GetLayer(0);

    // * * * * * * * * * * * * * * * * * * create Output Layer * * * * * * * * * * * * * * * * * *
    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ){
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }

    GDALDataset *poDS;
    poDS = poDriver->Create( tmpPath_clip.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ){
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer( "point_out", poMaskLayer->GetSpatialRef(), wkbLineString, NULL );
    if( poLayer == NULL ){
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    poDBLayer->Clip(poMaskLayer, poLayer);

    GDALClose( poDS );
    GDALClose( poDBLayer );
    GDALClose( poMaskLayer );
}

void bufferClip(const std::string& tmpPath_clip, const std::string& tmpPath_buffer){

    GDALDataset * poDS1;
    poDS1 = (GDALDataset*) GDALOpenEx( tmpPath_clip.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );

    if( poDS1 == NULL ){
        printf( "Open of clipped file failed.\n" );
        exit( 1 );
    }
    OGRLayer * poDBLayer;
    poDBLayer = poDS1->GetLayer(0);


    // * * * * * * * * * * * * * * * * * * create Output Layer * * * * * * * * * * * * * * * * * *
    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ){
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }

    GDALDataset *poDS;
    poDS = poDriver->Create( tmpPath_buffer.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ){
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer( "point_out", poDBLayer->GetSpatialRef(), wkbLineString, NULL );
    if( poLayer == NULL ){
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    OGRFieldDefn oField( "Angle", OFTReal );
    oField.SetPrecision(15);
    oField.SetWidth(32);

    if( poLayer->CreateField( &oField ) != OGRERR_NONE ){
        printf( "Creating Angle field failed.\n" );
        exit( 1 );
    }

    OGRFeature *poFeature;
    poDBLayer->ResetReading();

    // * * * * * * * * * * * * * * * * print out db content * * * * * * * * * * * * * * * *
    while( (poFeature = poDBLayer->GetNextFeature()) != NULL )
    {

        OGRFeatureDefn * poFDefn = poDBLayer->GetLayerDefn();

        OGRGeometry * poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        OGRLineString * poLineString;
        poLineString = (OGRLineString *) poGeometry;
        double Angle = 0.;
        for(unsigned int i = 0; i < poLineString->getNumPoints() - 1; i++){
            OGRPoint MyPoint1, MyPoint2;
            poLineString->getPoint(i, &MyPoint1);
            poLineString->getPoint(i + 1, &MyPoint2);
            computeAngle(Angle, MyPoint1, MyPoint2);

            OGRLineString LineString2;
            OGRLineString * poLineString2;
            LineString2.addPoint(MyPoint1.getX(), MyPoint1.getY());
            LineString2.addPoint(MyPoint2.getX(), MyPoint2.getY());
            poLineString2 = &LineString2;


            OGRGeometry * poBufferGeometry;
            poBufferGeometry = (OGRGeometry *) poLineString2;
            poBufferGeometry = poBufferGeometry->Buffer(.00005, 30);
            OGRLinearRing * poBufferPolygon = (OGRLinearRing *) poBufferGeometry;
            OGRFeature *poFeature;
            poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
            poFeature->SetField( "Angle", Angle );
            poBufferPolygon->closeRings();
            poFeature->SetGeometry( poBufferPolygon );

            if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
            {
                printf( "Failed to create feature in shapefile.\n" );
                exit( 1 );
            }

            OGRFeature::DestroyFeature( poFeature );

        }

    }

    GDALClose( poDS );
}

void rasterizeShp(const std::string& tmpPath_buffer, const std::string& tmpMaskImage){

}

int main(int argc, char *argv[])
{
    OGRRegisterAll();
    GDALAllRegister();

    const std::string AbsPathImage =      "/home/peter/Desktop/data/img/ON1678_clipped.tif";
    const std::string AbsPathDB =         "/home/peter/Desktop/data/GermanyRoads.sqlite";
    const std::string tmpPath_Rectangle = "/home/peter/Desktop/data/vec/rectangle.shp";
    const std::string tmpPath_clip =      "/home/peter/Desktop/data/vec/clip.shp";
    const std::string tmpPath_buffer =    "/home/peter/Desktop/data/vec/buffer.shp";

/*
    OGRLinearRing MyRing;
    OGRSpatialReference SRS;
    SRS.SetWellKnownGeogCS("WGS84");
    MyRing.assignSpatialReference(&SRS);

    createBoundingBox(AbsPathImage, MyRing);

    OGRPolygon MyPolygon;
    MyPolygon.assignSpatialReference(&SRS);
    MyPolygon.addRing(&MyRing);

    createTmpLayer(MyPolygon, tmpPath_Rectangle);

    clipDB(AbsPathDB, tmpPath_Rectangle, tmpPath_clip);
*/

    bufferClip(tmpPath_clip, tmpPath_buffer);


    //loadCompleteOSMToMem(MyPolygon);


    return 0;
}
