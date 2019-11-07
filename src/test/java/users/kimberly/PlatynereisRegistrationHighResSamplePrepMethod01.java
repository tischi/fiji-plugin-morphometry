package users.kimberly;

import de.embl.cba.morphometry.commands.PlatynereisRegistrationCommand;
import de.embl.cba.morphometry.registration.platynereis.PlatynereisRegistrationSettings;
import net.imagej.ImageJ;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class PlatynereisRegistrationHighResSamplePrepMethod01
{
	public static < R extends RealType< R > & NativeType< R > > void main( String[] args )
	{

		final ImageJ imageJ = new ImageJ();
		imageJ.ui().showUI();

		final PlatynereisRegistrationCommand< R > command = new PlatynereisRegistrationCommand<>();
		final PlatynereisRegistrationSettings settings = new PlatynereisRegistrationSettings();

		command.showIntermediateResults = true;
		command.outputResolution = 6;
		command.registrationResolution = settings.registrationResolution;
		command.opService = imageJ.op();
		command.inputDirectory = new File("/Users/tischer/Desktop/tomo_w7");
		command.outputDirectory = new File( "/Volumes/cba/exchange/Kimberly/High_res/samplePrepMethod01/tomo_w7-aligned" );
		command.inputResolutionMicrometer = 0.325;
		command.fileNameEndsWith = ".*.tiff";
		command.invertImage = true;

		command.run();

	}
}
