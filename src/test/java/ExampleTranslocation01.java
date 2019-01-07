import bdv.util.Bdv;
import bdv.util.BdvFunctions;
import bdv.util.BdvOptions;
import bdv.util.BdvStackSource;
import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Plots;
import de.embl.cba.morphometry.Rois;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.translocation.TranslocationComputer;
import de.embl.cba.morphometry.translocation.TranslocationResult;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.ARGBType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

import java.util.ArrayList;

public class ExampleTranslocation01
{



	public static  < T extends RealType< T > & NativeType< T > >
	void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		final ImagePlus imagePlus = IJ.openImage(
				ExampleTranslocation01.class.getResource(
						"translocation/test01.zip" ).getFile() );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", ExampleTranslocation01.class.getResource( "translocation/test01-rois.zip" ).getFile() );
		final ArrayList< FinalInterval > intervals = Rois.asIntervals( rm.getRoisAsArray() );

		ArrayList< RandomAccessibleInterval< T > > intensities = Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final TranslocationComputer computer = new TranslocationComputer(
				intensities,
				intervals,
				ij.op() );

		final ArrayList< TranslocationResult > results = computer.getResults();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks = createLabelMasks( intensities, results );

		final BdvStackSource< T > labelMasksSource = BdvFunctions.show(
				Views.stack( labelMasks ),
				"labelMasks",
				BdvOptions.options().is2D() );

		labelMasksSource.setDisplayRange( 0, results.size() + 2 );
		labelMasksSource.setColor( new ARGBType( ARGBType.rgba( 0,255,0,255 ) ) );

		final BdvStackSource< T > intensitiesSource = BdvFunctions.show(
				Views.stack( intensities ),
				"intensities",
				BdvOptions.options().addTo( labelMasksSource.getBdvHandle() ).is2D() );
		intensitiesSource.setDisplayRange( 0, Algorithms.getMaximumValue( intensities.get( 0 ) ) );


		for ( int r = 0; r < results.size(); r++ )
		{
			final CoordinatesAndValues coordinatesAndValues = new CoordinatesAndValues();

			for ( int t = 0; t < intensities.size(); t++ )
			{
				coordinatesAndValues.coordinates.add( 1.0 * t );
				coordinatesAndValues.values.add( ( double ) results.get( r ).translocation.get( t ) );
			}

			Plots.plot( coordinatesAndValues.coordinates, coordinatesAndValues.values, "time", "translocation - " + r );
		}


//		Utils.listOf2DImagesAsImagePlusMovie( labelMasks, "labelMasks" ).show();
//
//		for ( int i = 0; i < results.size(); i++ )
//		{
//			Utils.listOf2DImagesAsImagePlusMovie(
//					results.get( i ).cellMasks,
//					"cell masks 0" ).show();
//		}
//
//		imagePlus.show();

	}

	public static < T extends RealType< T > & NativeType< T > >
	ArrayList< RandomAccessibleInterval< T > > createLabelMasks(
			ArrayList< RandomAccessibleInterval< T > > movie,
			ArrayList< TranslocationResult > results )
	{
		final ArrayList< RandomAccessibleInterval< T > > labelMasks = new ArrayList<>();

		for ( int t = 0; t < movie.size(); t++ )
		{
			final RandomAccessibleInterval< T > labelMask = Utils.createEmptyArrayImg( movie.get( 0 ) );

			for ( int r = 0; r < results.size(); r++ )
			{
				RandomAccessibleInterval< T > mask = (RandomAccessibleInterval) results.get( r ).cellMasks.get( t );
				Utils.drawMaskIntoImage( mask, labelMask, r + 1 );
			}

			labelMasks.add( labelMask );

		}
		return labelMasks;
	}

}
