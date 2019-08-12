package tests;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.registration.drospholia.DrosophilaRegistrationSettings;
import de.embl.cba.morphometry.registration.drospholia.DrosophilaSingleChannelRegistration;
import ij.ImagePlus;
import loci.common.DebugTools;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.junit.Test;

import static de.embl.cba.morphometry.Utils.openWithBioFormats;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class TestDrosophilaRegistration < T extends RealType< T > & NativeType< T > >
{
	@Test
	public void testRegistration000()
	{

		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final OpService opService = ij.op();

		final String filePath = "../test-data/drosophila/low_res_x60_y55_z41_yaw-22.zip";
		final double[] actualCentre = { 60.0, 55.0, 41.0 };
		final int actualAngle = -22;

		final String inputPath =
				TestDrosophilaRegistration.class.getResource(
						filePath )
						.getFile().toString();

		final ImagePlus imagePlus = openWithBioFormats( inputPath );
		final double[] calibration = Utils.getCalibration( imagePlus );

		RandomAccessibleInterval< T > images =
				Utils.getChannelImages( imagePlus );

		runDrosophilaRegistrationTest( opService, actualCentre, actualAngle, calibration, images );
	}

	@Test
	public void testRegistration001()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final OpService opService = ij.op();

		final String filePath = "../test-data/drosophila/low_res_x58_y69_z38_yaw-54.zip";
		final double[] actualCentre = { 58.0, 69.0, 38.0 };
		final int actualAngle = -54;

		final String inputPath =
				TestDrosophilaRegistration.class.getResource(
						filePath )
						.getFile();

		final ImagePlus imagePlus = openWithBioFormats( inputPath );
		final double[] calibration = Utils.getCalibration( imagePlus );

		RandomAccessibleInterval< T > images =
				Utils.getChannelImages( imagePlus );

		runDrosophilaRegistrationTest(
				opService, actualCentre, actualAngle, calibration, images );
	}


	public void runDrosophilaRegistrationTest(
			OpService opService,
			double[] actualCentre,
			int actualAngle,
			double[] calibration,
			RandomAccessibleInterval< T > images )
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		RandomAccessibleInterval< T > image =
				Utils.getChannelImage( images, 0  );

		final DrosophilaRegistrationSettings settings =
				new DrosophilaRegistrationSettings();

		settings.onlyComputeEllipsoidParameters = true;

		final DrosophilaSingleChannelRegistration registration =
				new DrosophilaSingleChannelRegistration( settings, opService );

		registration.run( image, calibration );

		final double[] centre = registration.getElliposidCentreInInputImagePixelUnits();
		final double[] angles = registration.getElliposidEulerAnglesInDegrees();

		logCentre( actualCentre );
		logCentre( centre );
		logAngle( actualAngle );
		logAngle( angles[ 0 ] );

		assertArrayEquals( centre, actualCentre, 5.0 );
		assertEquals( angles[0], actualAngle, 5.0 );
	}

	public void logAngle( double angle )
	{
		Logger.log( "Angle: " + angle );
	}

	public void logCentre( double[] centre )
	{
		Logger.log( "Centre: " + centre[ 0 ] + ", " + centre[ 1 ] + ", " + centre[ 2 ] );
	}

}
