package explore;

import de.embl.cba.morphometry.ImageIO;
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

		String imageDataFolder = "/Users/tischer/Desktop/Hamburg_020_3_2_2k_Rec";
		String imageDataExtension = ".png";

		final ImagePlus imagePlus = FolderOpener.open( imageDataFolder, " file=(.*." + imageDataExtension + ")" );
		final double[] calibration = new double[]{1.0, 1.0, 1.0}; // TODO: get this information from somewhere

		final RandomAccessibleInterval< R > channelImages = ImageIO.getChannelImages( imagePlus );
		RandomAccessibleInterval< R > image = ImageIO.getChannelImage( channelImages, 0 );

		final PlatynereisRegistrationSettings settings = new PlatynereisRegistrationSettings();
		settings.showIntermediateResults = true;
		settings.registrationResolution = 8;

		final PlatynereisRegistration< R > registration = new PlatynereisRegistration<>( settings, opService );

		registration.run( image, calibration );
	}
}
