package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.registration.platynereis.PlatynereisRegistration;
import de.embl.cba.morphometry.registration.platynereis.PlatynereisRegistrationSettings;
import de.embl.cba.transforms.utils.Transforms;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.FolderOpener;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;
import java.util.ArrayList;


@Plugin(type = Command.class, menuPath = "Plugins>Registration>EMBL>Platynereis Registration" )
public class PlatynereisRegistrationCommand< R extends RealType< R > & NativeType< R > > implements Command
{
	@Parameter
	public UIService uiService;

	@Parameter
	public DatasetService datasetService;

	@Parameter
	public LogService logService;

	@Parameter
	public OpService opService;

	@Parameter
	public StatusService statusService;

	PlatynereisRegistrationSettings settings = new PlatynereisRegistrationSettings();

	@Parameter( style = "directory" )
	public File inputDirectory;

	@Parameter( style = "directory" )
	public File outputDirectory;

	@Parameter
	public String fileNameEndsWith = ".png";

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;

	@Parameter ( label = "Input image isotropic pixel size [micrometer]" )
	public double inputResolution;

	@Parameter
	public double outputResolution = settings.outputResolution;

	@Parameter
	public double registrationResolution = settings.registrationResolution;


	public void run()
	{
		setSettings();

		final double[] calibration = new double[]{
				inputResolution,
				inputResolution,
				inputResolution };

		// Open image
		//
		final ImagePlus imagePlus = FolderOpener.open( inputDirectory.getAbsolutePath(),
				" file=(.*" + fileNameEndsWith + ")" );
		final RandomAccessibleInterval< R > channelImages = ImageIO.getChannelImages( imagePlus );
		RandomAccessibleInterval< R > image = ImageIO.getChannelImage( channelImages, 0 );

		// Find registration
		//
		final PlatynereisRegistration< R > registration = new PlatynereisRegistration<>( settings, opService );
		registration.run( image, calibration );

		// Apply registration
		//
		Logger.log( "Creating registered image..." );
		ArrayList< RandomAccessibleInterval< R > > registeredImages =
				Transforms.transformAllChannels(
						channelImages,
						registration.getRegistrationTransform( calibration, outputResolution ),
						settings.getOutputImageInterval()
				);


		// Save registered image
		//
		final String outputFilePathStump = outputDirectory + File.separator + inputDirectory.getName();

		saveRegisteredImage( outputFilePathStump, registeredImages.get( 0 ) );

		Logger.log( "Done!" );
	}

	public void saveRegisteredImage( String outputFilePathStump, RandomAccessibleInterval< R > registeredImage )
	{

		// make 5D XYCZT
		registeredImage =
				Views.permute(
						Views.addDimension(
								Views.addDimension( registeredImage, 0 , 0 ),
								0 , 0 ),
						2 ,3 );

		final ImagePlus registered = ImageJFunctions.wrap( registeredImage, "transformed" );
		registered.getCalibration().setUnit( "micrometer" );
		registered.getCalibration().pixelWidth = settings.outputResolution;
		registered.getCalibration().pixelHeight = settings.outputResolution;
		registered.getCalibration().pixelDepth = settings.outputResolution;

		final String outputPath = outputFilePathStump + "-aligned.tif";
		Logger.log( "Saving registered image: " + outputPath );
		new FileSaver( registered ).saveAsTiff( outputPath );
	}


	public void setSettings()
	{
		settings.showIntermediateResults = showIntermediateResults;
		settings.registrationResolution = registrationResolution;
		settings.outputResolution = outputResolution;
	}


}
