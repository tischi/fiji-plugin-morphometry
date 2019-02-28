import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.drosophila.registration.DrosophilaRegistrationSettings;
import de.embl.cba.morphometry.drosophila.registration.DrosophilaSingleChannelRegistration;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidsMLJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import static de.embl.cba.morphometry.ImageIO.openWithBioFormats;

public class ExampleDrosophilaEllipsoidParameterComputation
{
	public static < T extends RealType< T > & NativeType< T > >
	void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final OpService opService = ij.op();

		final String inputPath =
				ExampleDrosophilaEllipsoidParameterComputation.class.getResource(
						"drosophila/test01_TR1_1_W0001_P0002_T0001.zip" ).getFile().toString();

		final ImagePlus imagePlus = openWithBioFormats( inputPath );

		final DrosophilaRegistrationSettings settings =
				new DrosophilaRegistrationSettings();
		settings.onlyComputeEllipsoidParameters = true;

		RandomAccessibleInterval< T > images =
				ImageIO.getChannelImages( imagePlus );

		RandomAccessibleInterval< T > image =
				ImageIO.getChannelImage( images, 0  );

		final DrosophilaSingleChannelRegistration registration =
				new DrosophilaSingleChannelRegistration( settings, opService );

		final double[] calibration = Utils.getCalibration( imagePlus );

		registration.run( image, calibration );

		final EllipsoidMLJ ellipsoidParameters = registration.getEllipsoidParameters();

		Logger.log( ellipsoidParameters.toString() );

		final AffineTransform3D alignmentTransform =
				EllipsoidsMLJ.createAlignmentTransform( ellipsoidParameters );

	}


}
