package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.fccf.FCCF;
import de.embl.cba.tables.Tables;
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

import javax.swing.*;
import java.io.File;
import java.util.*;

import static de.embl.cba.morphometry.fccf.FCCF.checkFile;

@Plugin( type = Command.class )
public class BDImageViewingAndSavingCommand extends DynamicCommand implements Initializable
{
	public static final String QC = "QC";
	public static final String PATH_PROCESSED_JPEG = "path_processed_jpeg";
	public static final String IMAGES_PROCESSED_JPEG = "images_processed_jpeg";
	public static File inputImagesDirectory;
	public static JTable jTable;
	public static File tableFile;
	public static String imagePathColumnName;
	public static String gateColumnName;

	private String[] fileNames;
	private int numFiles;

	@Parameter
	public LogService logService;

	@Parameter ( label = "Minimum File Size [kb]")
	public double minimumFileSizeKiloBytes = 10;

	@Parameter ( label = "Maximum File Size [kb]")
	public double maximumFileSizeKiloBytes = 100000;

	@Parameter ( label = "Minimum Brightfield Intensity" )
	public double minBF = 0.0;

	@Parameter ( label = "Maximum Brightfield Intensity" )
	public double maxBF = 1.0;

	@Parameter ( label = "Minimum GFP Intensity" )
	public double minGFP = 0.08;

	@Parameter ( label = "Maximum GFP Intensity" )
	public double maxGFP = 1.0;

	@Parameter ( label = "Processing Modality", choices = { FCCF.VIEW_RAW, FCCF.VIEW_PROCESSED_MONTAGE, FCCF.VIEW_PROCESSED_BF_GF_OVERLAY } )
	public String viewingModality = FCCF.VIEW_PROCESSED_MONTAGE;

	@Parameter ( label = "Preview Images from Gate", choices = {} )
	public String gateChoice = "";

	@Parameter ( label = "Preview Random Image", callback = "showRandomImage" )
	private Button showRandomImage;

	@Parameter ( label = "Output Directory" , style = "directory", persist = false)
	public File outputDirectory;

	@Parameter ( label = "Maximum Number of Files to Process" )
	private int maxNumFiles;

	private HashMap< String, ArrayList< Integer > > gateToRows;
	private ImagePlus rawImp;
	private ImagePlus processedImp;
	private int gateColumnIndex;
	private int pathColumnIndex;
	private static String experimentDirectory;
	private static String relativeProcessedImageDirectory;

	public void run()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		IJ.log( "\nSaving processed images: " );

		if ( jTable != null )
		{
			saveProcessedImagesAndTableWithQC();
		}
		else if ( inputImagesDirectory != null )
		{
			saveProcessedImages();
		}
	}

	public void saveProcessedImages()
	{
		final long currentTimeMillis = System.currentTimeMillis();

		for ( int i = 0; i < maxNumFiles; i++ )
		{
			final String inputImagePath = getInputImagePath( i );
			if ( checkFile( inputImagePath, minimumFileSizeKiloBytes, maximumFileSizeKiloBytes ) )
			{
				processedImp = createProcessedImagePlus( inputImagePath );
				saveImageAsJpeg( new File( inputImagePath ).getName(), processedImp );
			}
			Logger.progress( maxNumFiles, i + 1, currentTimeMillis, "Files saved" );
		}
	}

	public void saveProcessedImagesAndTableWithQC()
	{
		Tables.addColumn( jTable, QC, "Passed" );
		final int columnIndexQC = jTable.getColumnModel().getColumnIndex( QC );

		Tables.addColumn( jTable, PATH_PROCESSED_JPEG, "None" );
		final int columnIndexPath = jTable.getColumnModel().getColumnIndex( PATH_PROCESSED_JPEG );

		final long currentTimeMillis = System.currentTimeMillis();

		int rowCount = jTable.getRowCount();
		for ( int rowIndex = 0; rowIndex < rowCount; rowIndex++ )
		{
			final String inputImagePath = getInputImagePath( jTable, rowIndex );

			if ( ! checkFile( inputImagePath, minimumFileSizeKiloBytes, maximumFileSizeKiloBytes ) )
			{
				jTable.setValueAt( "Failed_FileSizeTooSmall", rowIndex, columnIndexQC );
				continue;
			}

			processedImp = createProcessedImagePlus( inputImagePath );

			final File path = saveImageAsJpeg( new File( inputImagePath ).getName(), processedImp );

			jTable.setValueAt( relativeProcessedImageDirectory + File.separator + path.getName(), rowIndex, columnIndexPath );

			Logger.progress( rowCount, rowIndex + 1, currentTimeMillis, "Files saved" );
		}

		saveTableWithAdditionalColumns();
	}

	public void saveTableWithAdditionalColumns()
	{
		IJ.log( "\nSaving table with additional columns..." );
		final File tableOutputFile = new File( tableFile.getAbsolutePath().replace( ".csv", "-processed-images.csv" ) );
		Tables.saveTable( jTable, tableOutputFile );
		IJ.log( "...done: " + tableOutputFile );
		IJ.log( " " );
		BDOpenTableCommand.glimpseTable( jTable );
	}

	@Override
	public void initialize()
	{
		getInfo(); // HACK: Workaround for bug in SJC.

		if ( jTable != null )
		{
			setGates();
			setProcessedJpegImagesOutputDirectory();
			pathColumnIndex = jTable.getColumnModel().getColumnIndex( imagePathColumnName );
		}
		else
		{
			setNoGateChoice();
			fetchFiles();
			maxNumFiles = fileNames.length;
		}
	}

	public void setProcessedJpegImagesOutputDirectory()
	{
		final MutableModuleItem<File> mutableInput = //
				getInfo().getMutableInput("outputDirectory", File.class);

		experimentDirectory = new File( tableFile.getParent() ).getParent();
		relativeProcessedImageDirectory = "../" + IMAGES_PROCESSED_JPEG;
		final File defaultValue = new File( experimentDirectory + File.separator +
				IMAGES_PROCESSED_JPEG );
		mutableInput.setValue( this, defaultValue );
		mutableInput.setDefaultValue( defaultValue );
	}

	public void setGates()
	{
		gateColumnIndex = jTable.getColumnModel().getColumnIndex( gateColumnName );
		gateToRows = Tables.uniqueColumnEntries( jTable, gateColumnIndex );

		final MutableModuleItem<String> gateChoiceItem = //
				getInfo().getMutableInput("gateChoice", String.class);

		gateChoiceItem.setChoices( new ArrayList<>( gateToRows.keySet() ) );
	}

	public void setNoGateChoice()
	{
		final MutableModuleItem<String> gateChoiceItem = //
				getInfo().getMutableInput("gateChoice", String.class);

		final ArrayList< String> choices = new ArrayList<>( );
		choices.add( "Option not available" );
		gateChoiceItem.setChoices( choices );
	}

	public void setMaxNumFiles()
	{
		final MutableModuleItem<Integer> item = //
				getInfo().getMutableInput("maxNumFiles", Integer.class);
		item.setDefaultValue( fileNames.length );
	}

	private void showRandomImage()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		if ( processedImp != null ) processedImp.close();

		final String filePath = getRandomFilePath();

		if ( ! checkFile( filePath, minimumFileSizeKiloBytes, maximumFileSizeKiloBytes ) ) return;

		processedImp = createProcessedImagePlus( filePath );

		processedImp.show();
	}

	private ImagePlus createProcessedImagePlus( String filePath )
	{
		ImagePlus processedImp = FCCF.createProcessedImage( filePath, getNameToRange(), FCCF.getNameToSlice(), viewingModality );
		return processedImp;
	}

	private HashMap< String, double[] > getNameToRange()
	{
		final HashMap< String, double[] > nameToRange = new HashMap<>();
		nameToRange.put( FCCF.BRIGHTFIELD, new double[]{ minBF, maxBF } );
		nameToRange.put( FCCF.GREEN_FLUORESCENCE, new double[]{ minGFP, maxGFP } );
		nameToRange.put( FCCF.FOREWARD_SCATTER, new double[]{ 0.0, 1.0 } );
		nameToRange.put( FCCF.SIDE_SCATTER, new double[]{ 0.0, 1.0 } );
		return nameToRange;
	}

	private String getRandomFilePath()
	{
		if ( jTable != null )
		{
			final Random random = new Random();
			final ArrayList< Integer > rowIndices = gateToRows.get( gateChoice );
			final Integer rowIndex = rowIndices.get( random.nextInt( rowIndices.size() ) );
			final String absoluteImagePath = getInputImagePath( jTable, rowIndex );
			return absoluteImagePath;
		}
		else if ( inputImagesDirectory != null )
		{
			fetchFiles();
			final Random random = new Random();
			final int randomInteger = random.nextInt( fileNames.length );
			final String absoluteImagePath = getInputImagePath( randomInteger );
			return absoluteImagePath;
		}
		else
		{
			return null;
		}
	}

	private String getInputImagePath( int randomInteger )
	{
		return inputImagesDirectory + File.separator + fileNames[ randomInteger ] ;
	}


	private void fetchFiles()
	{
		if ( fileNames != null ) return;
		fileNames = FCCF.readFileNamesFromDirectoryWithLogging( inputImagesDirectory );
	}

	private String getInputImagePath( JTable jTable, Integer rowIndex )
	{
		final String relativeImagePath = ( String ) jTable.getValueAt( rowIndex, pathColumnIndex );
		return tableFile.getParent() + File.separator + relativeImagePath;
	}



	private File saveImageAsJpeg( String fileName, ImagePlus outputImp )
	{
		final String outputFileName = fileName.replace( ".tiff", "" );
		final String outputFilePath = outputDirectory + File.separator + outputFileName + ".jpg";
		new File( outputDirectory.toString() ).mkdirs();
		new FileSaver( outputImp ).saveAsJpeg( outputFilePath  );
		return new File( outputFilePath );
	}


}
