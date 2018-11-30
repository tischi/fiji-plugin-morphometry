package de.embl.cba.morphometry.geometry.ellipsoids;

import ij.ImagePlus;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageByte;

public abstract class Ellipsoids3DImageSuite
{
	public static void fitEllipsoid( ImagePlus mask )
	{
		final ImageByte imageByte = new ImageByte( mask );
		Objects3DPopulation objects3DPopulation = new Objects3DPopulation( imageByte, 0 );
		final Object3D object = objects3DPopulation.getObject( 0 );
		System.out.println( "\n3D suite vectors");
		//System.out.println( object.getVolumePixels() );
		//System.out.println( object.getCenterAsVector() );
		System.out.println( object.getVectorAxis( 0 ) );
		System.out.println( object.getVectorAxis( 1 ) );
		System.out.println( object.getVectorAxis( 2 ) );
	}
}
