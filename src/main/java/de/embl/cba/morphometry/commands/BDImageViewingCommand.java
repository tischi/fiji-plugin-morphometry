package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.fccf.FCCF;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import loci.common.DebugTools;
import org.scijava.Initializable;
import org.scijava.command.Command;
import org.scijava.command.DynamicCommand;
import org.scijava.log.LogService;
import org.scijava.module.MutableModuleItem;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.widget.Button;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>FCCF>BD Processing" )
public class BDImageViewingCommand extends DynamicCommand implements Initializable
{
	@Parameter
	public LogService logService;

//	@Parameter ( style = "directory", label = "Output Directory" )
//	public File outputImageDirectory;

	@Parameter ( label = "Minimum File Size [kb]")
	public double minimumFileSizeKiloBytes = 10;

	@Parameter ( label = "Minimum Brightfield Intensity" )
	public double minBF = 0.08;

	@Parameter ( label = "Maximum Brightfield Intensity" )
	public double maxBF = 0.5;

	@Parameter ( label = "Minimum GFP Intensity" )
	public double minGFP = 0.08;

	@Parameter ( label = "Maximum GFP Intensity" )
	public double maxGFP = 1.0;

	@Parameter ( label = "Simple Overlay" )
	public boolean isSimpleOverlay = false;

	public static ArrayList< String > classChoices;

	@Parameter ( label = "Only Show Images from Class", choices = {} )
	public String classChoice = "";

	@Parameter ( label = "Preview Random Image", callback = "showRandomImage" )
	private Button showRandomImage;

	private String[] fileNames;
	private int numFiles;
	private ImagePlus inputImp;
	private ImagePlus outputImp;

	private ImagePlus brightFieldImp;
	private ImagePlus sideScatterImp;
	private ImagePlus forwardScatterImp;
	private ImagePlus gfpImp;
	private ImagePlus brightFieldGfpImp;

	public void run()
	{
		final long startMillis = System.currentTimeMillis();

		DebugTools.setRootLevel("OFF"); // Bio-Formats
//
//		fetchFiles();
//
//		createAndSaveImages();
//
//		logFinished( startMillis );
	}

	@Override
	public void initialize() {
		getInfo(); // HACK: Workaround for bug in SJC.

		final MutableModuleItem<String> classChoiceItem = //
				getInfo().getMutableInput("classChoice", String.class);

		classChoiceItem.setChoices( classChoices );
	}


//	private void showRandomImage()
//	{
//		if ( outputImp != null ) outputImp.close();
//
//		fetchFiles();
//
//		final String filePath = getRandomFilePath();
//
//		final ImagePlus processedImage = FCCF.createProcessedImage( filePath, minimumFileSizeKiloBytes, getNameToRange(), FCCF.getNameToSlice(), isSimpleOverlay );
//
//		if ( processedImage != null )
//		{
//			outputImp = processedImage;
//			outputImp.show();
//		}
//
//	}

	private HashMap< String, double[] > getNameToRange()
	{
		final HashMap< String, double[] > nameToRange = new HashMap<>();
		nameToRange.put( FCCF.BRIGHTFIELD, new double[]{ minBF, maxBF } );
		nameToRange.put( FCCF.GREEN_FLUORESCENCE, new double[]{ minGFP, maxGFP } );
		nameToRange.put( FCCF.FOREWARD_SCATTER, new double[]{ 0.0, 1.0 } );
		nameToRange.put( FCCF.SIDE_SCATTER, new double[]{ 0.0, 1.0 } );
		return nameToRange;
	}

//	private String getRandomFilePath()
//	{
//		final Random random = new Random();
//		return inputDirectory + File.separator + fileNames[ random.nextInt( numFiles ) ];
//	}
//
//	private void fetchFiles()
//	{
//		if ( fileNames != null ) return;
//
//		final long startMillis = System.currentTimeMillis();
//		IJ.log( "Fetching file list. Please wait..." );
//		fileNames = FCCF.getValidFileNames( inputDirectory );
//		IJ.log( "Fetched file list in " + ( System.currentTimeMillis() - startMillis) + " ms; number of files: " + numFiles );
//
//		numFiles = fileNames.length;
//	}

	private void showImages( String filePath, ImagePlus inputImp, ImagePlus outputImp )
	{
		final String fileName = new File( filePath ).getName();
		inputImp.setTitle( fileName );
		inputImp.show();
		outputImp.setTitle( fileName + "-output" );
		outputImp.show();
	}

//	private void saveImage( String fileName, ImagePlus outputImp )
//	{
//		final String outputFileName = fileName.replace( ".tiff", "" );
//		final String outputFilePath = outputImageDirectory + File.separator + outputFileName + ".jpg";
//		new File( outputImageDirectory.toString() ).mkdirs();
//		new FileSaver( outputImp ).saveAsJpeg( outputFilePath  );
//	}
//
//	private void logFinished( double startMillis )
//	{
//		IJ.log( "Processed " + numFiles + " in " + ( System.currentTimeMillis() - startMillis ) + " ms." );
//		Logger.log( "Done!" );
//	}
}
