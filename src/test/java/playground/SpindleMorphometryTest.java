package playground;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import de.embl.cba.morphometry.spindle.*;

import java.io.File;

public class SpindleMorphometryTest <T extends RealType<T> & NativeType< T > >
{

	public void run()
	{
		ImageJ imagej = new ImageJ();
		imagej.ui().showUI();

//		final String dapiPath = SpindleMorphometryTest.class.getResource( "spindle/test002-C01.tif" ).getFile();
//		final String tubulinPath = SpindleMorphometryTest.class.getResource( "spindle/test002-C00.tif" ).getFile();

		final String dapiPath = "";
		final String tubulinPath = "";

		// TODO: how to import data?

		IJ.run("Image Sequence...", "open=/Volumes/almfscreen/kletter/20180904/Subfolder/experiment--2018_09_04_15_56_43/CAM1--2018_09_04_18_07_15/slide--S00/chamber--U00--V00/field--X00--Y00 file=(.*J28.*C00.*) sort");
		final ImagePlus tubulinImp = IJ.getImage();

		IJ.run("Image Sequence...", "open=/Volumes/almfscreen/kletter/20180904/Subfolder/experiment--2018_09_04_15_56_43/CAM1--2018_09_04_18_07_15/slide--S00/chamber--U00--V00/field--X00--Y00 file=(.*J28.*C01.*) sort");
		final ImagePlus dapiImp = IJ.getImage();

		String outputDirectory = "/Volumes/almfscreen/kletter/20180904/Subfolder-out";

		final double[] calibration = Utils.getCalibration( dapiImp );

		final Img< T > dapi = ImageJFunctions.wrapReal( dapiImp );
		final Img< T > tubulin = ImageJFunctions.wrapReal( tubulinImp );

		SpindleMorphometrySettings settings = new SpindleMorphometrySettings();
		settings.showIntermediateResults = true;
		settings.inputCalibration = Utils.getCalibration( dapiImp );
		settings.dnaImage = dapi;
		settings.tubulinImage = tubulin;
		settings.workingVoxelSize = 0.25;
		settings.maxDnaAxisDist = 6;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
		settings.watershedSeedsGlobalDistanceThreshold = 2.0;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = new File( outputDirectory );
		settings.inputDataSetName = "test";


		SpindleMorphometry morphometry = new SpindleMorphometry( settings, imagej.op() );
		morphometry.run();

	}

	public static void main( String... args )
	{
		new SpindleMorphometryTest().run();

	}
}
