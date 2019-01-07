package de.embl.cba.morphometry.translocation;

import bdv.util.BdvFunctions;
import bdv.util.BdvOptions;
import bdv.util.BdvStackSource;
import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.ARGBType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.display.DisplayService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.io.File;
import java.util.ArrayList;

@Plugin(type = Command.class, menuPath = "Plugins>Measure Plasma Membrane Translocation" )
public class TranslocationCommand<T extends RealType<T> & NativeType< T > > implements Command
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
		ArrayList< RandomAccessibleInterval< T > > intensities = getIntensitiesFromImagePlus();

		final ArrayList< FinalInterval > intervals = getIntervalsFromRoiManager();

		final TranslocationComputer computer =
				new TranslocationComputer(
					intensities,
					intervals,
					opService );

		final ArrayList< TranslocationResult > results = computer.getResults();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks = createLabelMasks( intensities, results );

		showLabelMasksAndIntensities( intensities, labelMasks );

		plotTranslocations( results );
	}

	public ArrayList< RandomAccessibleInterval< T > > getIntensitiesFromImagePlus()
	{
		return Utils.get2DImagePlusMovieAsFrameList(
				imagePlus,
				1 );
	}

	public ArrayList< FinalInterval > getIntervalsFromRoiManager()
	{
		final RoiManager rm = RoiManager.getInstance();

		if ( rm == null || rm.getRoisAsArray().length == 0 )
		{
			Logger.error( "This plugin requires rectangular ROIs added to the ROIManager." );
			return null;
		}

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

	public static < T extends RealType< T > & NativeType< T > > RandomAccessibleInterval< T > createLabelMask( ArrayList< RandomAccessibleInterval< T > > movie, ArrayList< TranslocationResult > results, int t )
	{
		final RandomAccessibleInterval< T > labelMask = Utils.createEmptyArrayImg( movie.get( 0 ) );

		for ( int r = 0; r < results.size(); r++ )
		{
			RandomAccessibleInterval< T > mask = (RandomAccessibleInterval) results.get( r ).cellMasks.get( t );
			Utils.drawMaskIntoImage( mask, labelMask, r + 1 );
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
					"translocation-" + ( region + 1 ) );
		}
	}

	public static CoordinatesAndValues getTranslocations( ArrayList< TranslocationResult > results, int r )
	{
		final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

		for ( int t = 0; t < results.get( r ).insideIntensities.size(); t++ )
		{
			coordinatesAndValues.coordinates.add( 1.0 * t );
			coordinatesAndValues.values.add( ( double ) results.get( r ).translocation.get( t ) );
		}
		return coordinatesAndValues;
	}

	public static < T extends RealType< T > & NativeType< T > >
	void showLabelMasksAndIntensities(
			ArrayList< RandomAccessibleInterval< T > > intensities,
			ArrayList< RandomAccessibleInterval< T > > labelMasks )
	{
		final BdvStackSource< T > labelMasksSource = BdvFunctions.show(
				Views.stack( labelMasks ),
				"labelMasks",
				BdvOptions.options().is2D() );

		labelMasksSource.setDisplayRange( 0, Algorithms.getMaximumValue( labelMasks.get( 0 ) ) + 2 );
		labelMasksSource.setColor( new ARGBType( ARGBType.rgba( 0,255,0,255 ) ) );

		final BdvStackSource< T > intensitiesSource = BdvFunctions.show(
				Views.stack( intensities ),
				"intensities",
				BdvOptions.options().addTo( labelMasksSource.getBdvHandle() ).is2D() );
		intensitiesSource.setDisplayRange( 0, Algorithms.getMaximumValue( intensities.get( 0 ) ) );
	}


}
