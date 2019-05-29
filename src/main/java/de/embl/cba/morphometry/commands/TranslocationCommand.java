package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Plots;
import de.embl.cba.morphometry.Rois;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.translocation.MembraneTranslocationComputer;
import de.embl.cba.morphometry.translocation.TranslocationResult;
import de.embl.cba.tables.ExploreIntensityImageAndLabelImageAndTable;
import de.embl.cba.tables.Tables;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.frame.RoiManager;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.util.ArrayList;

@Plugin(type = Command.class, menuPath = "Plugins>Measure>Membrane Translocation" )
public class TranslocationCommand< R extends RealType< R > & NativeType< R > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter
	public ImagePlus intensitiesImp;

	@Parameter( label = "Output directory",  style = "directory" )
	public File outputDirectory;

	@Parameter ( label = "Show intermediate results" )
	public boolean showIntermediateResults = false;

	@Parameter ( label = "Show translocation plots" )
	public boolean showTranslocationPlots = false;

	@Parameter ( label = "Review membrane segmentation" )
	public boolean reviewMembraneSegmentation = false;
	private String membraneLabelsPath;
	private String intensitiesPath;
	private String tablePath;
	private String outputPathStump;


	public void run()
	{
		// fetch input
		ArrayList< RandomAccessibleInterval< R > > intensities =
				Utils.get2DImagePlusMovieAsFrameList(
						intensitiesImp,
						1 );

		final ArrayList< FinalInterval > intervals =
				getIntervalsFromRoiManagerAndSaveRois();

		// compute
		final MembraneTranslocationComputer computer =
				new MembraneTranslocationComputer(
					intensities,
					intervals,
					opService );

		final ArrayList< TranslocationResult > results = computer.getResults();

		// process output

		outputPathStump = outputDirectory + File.separator + intensitiesImp.getTitle();
		intensitiesPath = outputPathStump + "-intensities.tif";
		membraneLabelsPath = outputPathStump + "-membrane-labels.tif";

		final ImagePlus labelMasksImp =
				Utils.getAsImagePlusMovie(
						createLabelMasks( intensities, results ),
						"label masks" );

		final JTable jTable = TranslocationResult.resultsAsTable( results );

		Tables.addColumn( jTable,
				"Path_Membrane_Labels",
				new File( membraneLabelsPath ).getName() );

		Tables.addColumn( jTable,
				"Path_Intensity_Image",
				new File( intensitiesPath ).getName() );

		saveResults( jTable, labelMasksImp, intensitiesImp );

		if ( showTranslocationPlots )
			plotTranslocations( results );

		if ( reviewMembraneSegmentation )
		{
			new ExploreIntensityImageAndLabelImageAndTable(
					intensitiesPath,
					membraneLabelsPath,
					tablePath,
					true,
					false
			);
		}


	}

	private void saveResults(
			JTable jTable,
			ImagePlus labelMasksImp,
			ImagePlus intensitiesImp )
	{
		new FileSaver( labelMasksImp ).saveAsTiff( membraneLabelsPath );
		new FileSaver( intensitiesImp ).saveAsTiff( intensitiesPath );
		tablePath = outputPathStump + "-table.csv";
		Tables.saveTable( jTable, new File( tablePath ) );
	}


	private ArrayList< FinalInterval > getIntervalsFromRoiManagerAndSaveRois()
	{
		final RoiManager rm = RoiManager.getInstance();

		if ( rm == null || rm.getRoisAsArray().length == 0 )
		{
			Logger.error( "This plugin requires rectangular ROIs added to the ROIManager." );
			return null;
		}

		// select all
		for ( int i = 0; i < rm.getCount(); i++ )
			rm.select( i );

		rm.runCommand( "save",
				outputDirectory
						+ File.separator
						+ intensitiesImp.getTitle()
						+ "-rois.zip"  );

		return Rois.asIntervals( rm.getRoisAsArray() );
	}

	private static < T extends RealType< T > & NativeType< T > >
	ArrayList< RandomAccessibleInterval< T > > createLabelMasks(
			ArrayList< RandomAccessibleInterval< T > > movie,
			ArrayList< TranslocationResult > results )
	{
		final ArrayList< RandomAccessibleInterval< T > > labelMasks = new ArrayList<>();

		for ( int t = 0; t < movie.size(); t++ )
		{
			final RandomAccessibleInterval< T > labelMask =
					createLabelMask( movie, results, t );
			labelMasks.add( labelMask );
		}

		return labelMasks;
	}

	private static < T extends RealType< T > & NativeType< T > > RandomAccessibleInterval< T >
	createLabelMask(
			ArrayList< RandomAccessibleInterval< T > > movie,
			ArrayList< TranslocationResult > results,
			int t )
	{
		final RandomAccessibleInterval< T > labelMask = Utils.createEmptyCopy( movie.get( 0 ) );

		for ( int region = 0; region < results.size(); region++ )
		{
			Utils.drawMaskIntoImage(
					(RandomAccessibleInterval) results.get( region ).membranes.get( t ),
					labelMask,
					region + 1 );

			for ( int i = 0; i < 2; i++ )
			{
				Utils.drawMaskIntoImage(
						( RandomAccessibleInterval ) ( ( ArrayList ) results.get( region ).nonMembraneMasks.get( t ) ).get( i ),
						labelMask,
						region + 1 );
			}

		}

		return labelMask;
	}

	private static void plotTranslocations( ArrayList< TranslocationResult > results )
	{
		for ( int region = 0; region < results.size(); region++ )
		{
			Plots.plot(
					getTranslocations( results, region ),
					"time",
					"translocation - region " + ( region + 1 ) );
		}
	}

	private static CoordinatesAndValues getTranslocations(
			ArrayList< TranslocationResult > results, int r )
	{
		final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

		for ( int t = 0; t < results.get( r ).brighterIntensities.size(); t++ )
		{
			if ( results.get( r ).translocations.get( t ) != null )
			{
				coordinatesAndValues.coordinates.add( 1.0 * t );
				coordinatesAndValues.values.add(
						( Double ) results.get( r ).translocations.get( t ) );
			}
		}
		return coordinatesAndValues;
	}

}
