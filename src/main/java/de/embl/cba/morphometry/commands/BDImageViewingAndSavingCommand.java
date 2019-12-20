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
	@Parameter
	public LogService logService;

//	@Parameter ( style = "directory", label = "Output Directory" )
//	public File outputImageDirectory;

	@Parameter ( label = "Minimum File Size [kb]")
	public double minimumFileSizeKiloBytes = 10;

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

	// TODO: Could this also be (resolved) parameters?
	public static JTable jTable;
	public static File tableFile;
	public static String imagePathColumnName;
	public static String gateColumnName;

	private HashMap< String, ArrayList< Integer > > gateToRows;
	private ImagePlus rawImp;
	private ImagePlus processedImp;
	private int gateColumnIndex;
	private int pathColumnIndex;

	public void run()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats
		saveProcessedImagesAndTableWithQC();
	}

	public void saveProcessedImagesAndTableWithQC()
	{

		IJ.log( "\nSaving processed images: " );
		Tables.addColumn( jTable, QC, "Passed" );
		final int columnIndexQC = jTable.getColumnModel().getColumnIndex( QC );

		Tables.addColumn( jTable, PATH_PROCESSED_JPEG, "None" );
		final int columnIndexPath = jTable.getColumnModel().getColumnIndex( PATH_PROCESSED_JPEG );

		final long currentTimeMillis = System.currentTimeMillis();

		int rowCount = jTable.getRowCount();
		for ( int rowIndex = 0; rowIndex < rowCount; rowIndex++ )
		{
			final String inputImagePath = getInputImagePath( jTable, rowIndex );

			if ( ! checkFile( inputImagePath, minimumFileSizeKiloBytes ) )
			{
				jTable.setValueAt( "Failed_FileSizeTooSmall", rowIndex, columnIndexQC );
				continue;
			}

			processedImp = createProcessedImagePlus( inputImagePath );

			final File path = saveImageAsJpeg( new File( inputImagePath ).getName(), processedImp );
			jTable.setValueAt( path.getAbsolutePath(), rowIndex, columnIndexPath );

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
		setGates();
		setProcessedJpegImagesOutputDirectory();
		pathColumnIndex = jTable.getColumnModel().getColumnIndex( imagePathColumnName );
	}

	public void setProcessedJpegImagesOutputDirectory()
	{
		final MutableModuleItem<File> mutableInput = //
				getInfo().getMutableInput("outputDirectory", File.class);

		final String experimentDirectory = new File( tableFile.getParent() ).getParent();
		final File defaultValue = new File( experimentDirectory + File.separator + "images_processed_jpeg" );
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

	private void showRandomImage()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		if ( processedImp != null ) processedImp.close();

		final String filePath = getRandomFilePath();

		if ( ! checkFile( filePath, minimumFileSizeKiloBytes ) ) return;

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
		final Random random = new Random();
		final ArrayList< Integer > rowIndices = gateToRows.get( gateChoice );
		final Integer rowIndex = rowIndices.get( random.nextInt( rowIndices.size() ) );

		final String absoluteImagePath = getInputImagePath( jTable, rowIndex );

		return absoluteImagePath;
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
