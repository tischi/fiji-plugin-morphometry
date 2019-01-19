package de.embl.cba.morphometry.commands;

import bdv.util.*;
import de.embl.cba.bdv.utils.BdvUtils;
import de.embl.cba.bdv.utils.overlays.BdvGrayValuesOverlay;
import de.embl.cba.bdv.utils.selection.BdvSelectionEventHandler;
import de.embl.cba.bdv.utils.sources.SelectableARGBConvertedRealSource;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.translocation.MembraneTranslocationComputer;
import de.embl.cba.morphometry.translocation.TranslocationResult;
import de.embl.cba.tables.TableUtils;
import de.embl.cba.tables.objects.ObjectTablePanel;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.frame.RoiManager;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.display.DisplayService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

@Plugin(type = Command.class, menuPath = "Plugins>Measurement>Measure Plasma Membrane Translocation" )
public class TranslocationCommand< R extends RealType< R > & NativeType< R > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter
	public DisplayService displayService;

	@Parameter
	public ImagePlus imagePlus;

	@Parameter( style = "directory" )
	public File outputDirectory;

	@Parameter
	public boolean showIntermediateResults = false;

	public void run()
	{
		// fetch input
		ArrayList< RandomAccessibleInterval< R > > intensities = getIntensitiesFromImagePlus();
		final ArrayList< FinalInterval > intervals = getIntervalsFromRoiManagerAndSaveRois();
		String dataSetName = imagePlus.getTitle();

		final MembraneTranslocationComputer computer =
				new MembraneTranslocationComputer(
					intensities,
					intervals,
					opService );

		final ArrayList< TranslocationResult > results = computer.getResults();

		final ArrayList< RandomAccessibleInterval< R > > labelMasks = createLabelMasks( intensities, results );

		final Bdv bdv = showLabelMasksAndIntensities( intensities, labelMasks, dataSetName );

		plotTranslocations( results );

		final ObjectTablePanel objectTablePanel = showResultsAsTable( results );

		final ImagePlus labelMasksImp = Utils.listOf2DImagesAsImagePlusMovie( labelMasks, "label masks" );

		saveResults( objectTablePanel, labelMasksImp );

		// TODO: connect table and bdv

	}

	public void saveResults( ObjectTablePanel objectTablePanel, ImagePlus labelMasksImp )
	{
		final String outputPathStump = outputDirectory + File.separator + imagePlus.getTitle();

		new FileSaver( labelMasksImp ).saveAsTiff( outputPathStump  + "-labels.tif" );

		try
		{
			TableUtils.saveTable( objectTablePanel.getTable(), new File( outputPathStump + "-table.csv") );
		}
		catch ( IOException e )
		{
			e.printStackTrace();
		}
	}

	public ObjectTablePanel showResultsAsTable( ArrayList< TranslocationResult > results )
	{
		final JTable table = TranslocationResult.resultsAsTable( results );
		final ObjectTablePanel objectTablePanel = new ObjectTablePanel( table );
		objectTablePanel.showPanel();

		return objectTablePanel;
	}

	public ArrayList< RandomAccessibleInterval< R > > getIntensitiesFromImagePlus()
	{
		return Utils.get2DImagePlusMovieAsFrameList(
				imagePlus,
				1 );
	}

	public ArrayList< FinalInterval > getIntervalsFromRoiManagerAndSaveRois()
	{
		final RoiManager rm = RoiManager.getInstance();

		if ( rm == null || rm.getRoisAsArray().length == 0 )
		{
			Logger.error( "This plugin requires rectangular ROIs added to the ROIManager." );
			return null;
		}

		// select all
		for ( int i = 0; i < rm.getCount(); i++ )
		{
			rm.select( i );
		}
		rm.runCommand( "save", outputDirectory + File.separator + imagePlus.getTitle() + "-rois.zip"  );

		return Rois.asIntervals( rm.getRoisAsArray() );
	}

	public static < T extends RealType< T > & NativeType< T > >
	ArrayList< RandomAccessibleInterval< T > > createLabelMasks(
			ArrayList< RandomAccessibleInterval< T > > movie,
			ArrayList< TranslocationResult > results )
	{
		final ArrayList< RandomAccessibleInterval< T > > labelMasks = new ArrayList<>();

		for ( int t = 0; t < movie.size(); t++ )
		{
			final RandomAccessibleInterval< T > labelMask = createLabelMask( movie, results, t );

			labelMasks.add( labelMask );
		}

		return labelMasks;
	}

	public static < T extends RealType< T > & NativeType< T > > RandomAccessibleInterval< T >
	createLabelMask( ArrayList< RandomAccessibleInterval< T > > movie, ArrayList< TranslocationResult > results, int t )
	{
		final RandomAccessibleInterval< T > labelMask = Utils.createEmptyCopy( movie.get( 0 ) );

		for ( int region = 0; region < results.size(); region++ )
		{
			Utils.drawMaskIntoImage(
					(RandomAccessibleInterval) results.get( region ).membranes.get( t ),
					labelMask,
					region + 1 );

			Utils.drawMaskIntoImage(
					(RandomAccessibleInterval) results.get( region ).insideOutsideMasks.get( t ),
					labelMask,
					region + 1 );
		}

		return labelMask;
	}

	public static void plotTranslocations( ArrayList< TranslocationResult > results )
	{
		for ( int region = 0; region < results.size(); region++ )
		{
			Plots.plot(
					getTranslocations( results, region ),
					"time",
					"translocation - region " + ( region + 1 ) );
		}
	}

	public static CoordinatesAndValues getTranslocations( ArrayList< TranslocationResult > results, int r )
	{
		final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

		for ( int t = 0; t < results.get( r ).brighterIntensities.size(); t++ )
		{
			if ( results.get( r ).translocations.get( t ) != null )
			{
				coordinatesAndValues.coordinates.add( 1.0 * t );
				coordinatesAndValues.values.add( ( Double ) results.get( r ).translocations.get( t ) );
			}
		}
		return coordinatesAndValues;
	}

	public static < R extends RealType< R > & NativeType< R > >
	BdvHandle showLabelMasksAndIntensities(
			ArrayList< RandomAccessibleInterval< R > > intensities,
			ArrayList< RandomAccessibleInterval< R > > labelMasks,
			String dataSetName )
	{
		final RandomAccessibleIntervalSource4D< R > labelsSource =
				BdvUtils.createSourceFrom2DFrameList( labelMasks, "labels" );

		final SelectableARGBConvertedRealSource argbLabelsSource =
				new SelectableARGBConvertedRealSource( labelsSource );

		final BdvStackSource< R > labelMasksSource = BdvFunctions.show(
				argbLabelsSource,
				labelMasks.size(),
				BdvOptions.options().is2D() );

		final BdvHandle bdvHandle = labelMasksSource.getBdvHandle();

		new BdvGrayValuesOverlay( bdvHandle, 25);

		new BdvSelectionEventHandler( bdvHandle, argbLabelsSource );

		final BdvStackSource< R > intensitiesSource = BdvFunctions.show(
				Views.stack( intensities ),
				"intensities",
				BdvOptions.options().addTo( bdvHandle ).is2D() );

		intensitiesSource.setDisplayRange( 0, Algorithms.getMaximumValue( intensities.get( 0 ) ) );

		return bdvHandle;
	}


}
