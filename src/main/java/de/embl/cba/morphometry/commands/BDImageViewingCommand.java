package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.fccf.FCCF;
import de.embl.cba.tables.Tables;
import ij.ImagePlus;
import loci.common.DebugTools;
import org.scijava.Initializable;
import org.scijava.command.Command;
import org.scijava.command.DynamicCommand;
import org.scijava.command.Interactive;
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
public class BDImageViewingCommand extends DynamicCommand implements Initializable, Interactive
{
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

	@Parameter ( label = "Viewing Modality", choices = { FCCF.VIEW_RAW, FCCF.VIEW_PROCESSED_MONTAGE, FCCF.VIEW_PROCESSED_BF_GF_OVERLAY } )
	public String viewingModality = FCCF.VIEW_PROCESSED_MONTAGE;

	@Parameter ( label = "Only Show Images from Gate", choices = {} )
	public String gateChoice = "";

	@Parameter ( label = "Preview Random Image", callback = "showRandomImage" )
	private Button showRandomImage;

	// TODO: Could this also be (resolved) parameters?
	public static JTable jTable;
	public static String imagePathColumnName;
	public static String gateColumnName;
	public static String imagesRootDir;

	private HashMap< String, ArrayList< Integer > > gateToRows;
	private ImagePlus rawImp;
	private ImagePlus processedImp;
	private int gateColumnIndex;
	private int pathColumnIndex;

	public void run()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats
//
//		fetchFiles();
//
//		createAndSaveImages();
//
//		logFinished( startMillis );
	}

	@Override
	public void initialize()
	{
		getInfo(); // HACK: Workaround for bug in SJC.
		setGates();
		pathColumnIndex = jTable.getColumnModel().getColumnIndex( imagePathColumnName );
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
		if ( processedImp != null ) processedImp.close();

		final String filePath = getRandomFilePath();

		if ( ! checkFile( filePath, minimumFileSizeKiloBytes ) ) return;

		processedImp = FCCF.createProcessedImage( filePath, getNameToRange(), FCCF.getNameToSlice(), viewingModality );

		processedImp.show();
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

		final String relativeImagePath = ( String ) jTable.getValueAt( rowIndex, pathColumnIndex );

		return imagesRootDir + File.separator + relativeImagePath;
	}
}
