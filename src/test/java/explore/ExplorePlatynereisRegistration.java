package explore;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.registration.platynereis.PlatynereisRegistration;
import de.embl.cba.morphometry.registration.platynereis.PlatynereisRegistrationSettings;
import ij.ImagePlus;
import ij.plugin.FolderOpener;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

public class ExplorePlatynereisRegistration
{
	public static < R extends RealType< R > & NativeType< R > > void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		final OpService opService = ij.op();

		String imageDataExtension = ".png";

//		String imageDataFolder = "/Volumes/cba/exchange/Kimberly/Low_res/020_3_2__w7";
//		final double[] calibration = new double[]{1.0, 1.0, 1.0};

//		String imageDataFolder = "/Volumes/cba/exchange/Kimberly/Low_res/038_2__w2";
//		final double[] calibration = new double[]{1.0, 1.0, 1.0};

		String imageDataFolder = "/Volumes/cba/exchange/Kimberly/Low_res/PLATY3A";
		final double[] calibration = new double[]{0.70022, 0.70022, 0.70022};

//		String imageDataFolder = "/Volumes/cba/exchange/Kimberly/Low_res/PLATY3B";
//		final double[] calibration = new double[]{0.70022, 0.70022, 0.70022};

		final ImagePlus imagePlus = FolderOpener.open( imageDataFolder, " file=(.*." + imageDataExtension + ")" );

		final RandomAccessibleInterval< R > channelImages = Utils.getChannelImages( imagePlus );
		RandomAccessibleInterval< R > image = Utils.getChannelImage( channelImages, 0 );

		final PlatynereisRegistrationSettings settings = new PlatynereisRegistrationSettings();
		settings.showIntermediateResults = true;
		settings.registrationResolution = 8; // micrometer
		settings.inputCalibration = calibration;

		final PlatynereisRegistration< R > registration = new PlatynereisRegistration<>( settings, opService );

		registration.run( image);
	}
}
