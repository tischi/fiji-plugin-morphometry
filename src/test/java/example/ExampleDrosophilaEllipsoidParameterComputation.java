package example;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.registration.drospholia.DrosophilaRegistrationSettings;
import de.embl.cba.morphometry.registration.drospholia.DrosophilaSingleChannelRegistration;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import static de.embl.cba.morphometry.Utils.openWithBioFormats;

public class ExampleDrosophilaEllipsoidParameterComputation
{
	public static < T extends RealType< T > & NativeType< T > >
	void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final OpService opService = ij.op();

//		final String inputPath =
//				ExampleDrosophilaEllipsoidParameterComputation.class.getResource(
//						"../drosophila/test01_TR1_1_W0001_P0002_T0001.zip" ).getFile().toString();

		final String inputPath =
				"/Users/tischer/Desktop/tim-wrong/test03_TR1_1_W0001_P0003_T0001.lsm";

		final ImagePlus imagePlus = openWithBioFormats( inputPath );

		imagePlus.show();

		final DrosophilaRegistrationSettings settings =
				new DrosophilaRegistrationSettings();
		settings.onlyComputeEllipsoidParameters = true;
//		settings.showIntermediateResults = true;

		RandomAccessibleInterval< T > images =
				Utils.getChannelImages( imagePlus );

		RandomAccessibleInterval< T > image =
				Utils.getChannelImage( images, 0  );

		final DrosophilaSingleChannelRegistration registration =
				new DrosophilaSingleChannelRegistration( settings, opService );

		final double[] calibration = Utils.getCalibration( imagePlus );

		if ( ! registration.run( image, calibration ) )
		{
			Logger.log( "Error: Could not segment embryo." );
		}
		else
		{
			final double[] centre = registration.getElliposidCentreInInputImagePixelUnits();
			final double[] angles = registration.getElliposidEulerAnglesInDegrees();
			Logger.log( "Centre: " + centre[ 0 ] + ", " + centre[ 1 ] + ", " + centre[ 2 ] );
			Logger.log( "Angles: " + angles[ 0 ] + ", " + angles[ 1 ] + ", " + angles[ 2 ] );
		}
	}


}
