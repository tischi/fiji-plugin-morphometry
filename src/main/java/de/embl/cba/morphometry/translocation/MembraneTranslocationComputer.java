package de.embl.cba.morphometry.translocation;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.regions.RegionAndSize;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.morphometry.segmentation.SignalOverBackgroundSegmenter;
import de.embl.cba.morphometry.skeleton.Skeletons;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.converter.Converters;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.view.Views;

import java.util.ArrayList;
import java.util.Collections;

import static net.imglib2.algorithm.labeling.ConnectedComponents.StructuringElement.FOUR_CONNECTED;

public class MembraneTranslocationComputer< R extends RealType< R > & NativeType< R > >
{
	public static final Double SEGMENTATION_ERROR = null;
	final ArrayList< RandomAccessibleInterval< R > > inputs;
	final ArrayList<FinalInterval> intervals;
	private double signalToNoise;
	private int closingRadius;
	private int resolutionBlurWidthInPixel;
	private int blurRadius;
	private int erosionRadius;
	private int region;

	public ArrayList< TranslocationResult > getResults()
	{
		return results;
	}

	final ArrayList< TranslocationResult > results;

	final OpService opService;
	private int minimalRegionSize;

	public MembraneTranslocationComputer(
			ArrayList< RandomAccessibleInterval< R > > inputs,
			ArrayList< FinalInterval > intervals,
			OpService opService )
	{
		this.inputs = inputs;
		this.intervals = intervals;
		this.opService = opService;

		minimalRegionSize = 50;
		signalToNoise = 25.0;
		closingRadius = 5;
		resolutionBlurWidthInPixel = 3;
		blurRadius = 4;
		erosionRadius = 4;

		results = new ArrayList<TranslocationResult >();

		computeTranslocations();
	}

	private void computeTranslocations()
	{
		for ( region = 0; region < intervals.size(); region++ )
			computeTranslocationsInRoi( intervals.get( region ) );

	}

	private void computeTranslocationsInRoi( FinalInterval interval )
	{
		final TranslocationResult< ? > result = new TranslocationResult();

		final ArrayList nonMembraneMasks = result.nonMembraneMasks;

		measureRegionCenter( interval, result );

		for ( int t = 0; t < inputs.size(); t++ )
		{
			createMembraneAndOtherRegions( interval, result, t );

			measureIntensities( result, t  );

			computeTranslocation( result, t );
		}

		results.add( result );
	}

	private void measureRegionCenter( FinalInterval interval, TranslocationResult result )
	{
		result.regionCenterX = Utils.getCenterLocation( interval )[ 0 ];
		result.regionCenterY = Utils.getCenterLocation( interval )[ 1 ];
	}

	private void computeTranslocation( TranslocationResult result, int t )
	{
		Double membrane = ( Double ) result.membraneIntensities.get( t );
		Double outside = ( Double ) result.dimmerIntensities.get( t );
		Double inside = ( Double ) result.brighterIntensities.get( t );

		if ( inside != null && outside != null && membrane != null )
		{
			final double translocation = ( membrane - outside ) / ( inside - outside );
			result.translocations.add( translocation );
		}
		else
		{
			// something went wrong
			result.translocations.add( SEGMENTATION_ERROR );
		}
	}

	private void createMembraneAndOtherRegions(
			FinalInterval interval,
			TranslocationResult result,
			int t )
	{
		final RandomAccessibleInterval< R > intensity =
				Views.interval( inputs.get( t ) , interval );

		RandomAccessibleInterval< BitType > cellEdgeSkeleton =
				computeCellEdge( result, intensity );

		final RandomAccessibleInterval< BitType > membrane =
				Skeletons.longestBranch( cellEdgeSkeleton );

		result.cellEdges.add( cellEdgeSkeleton );
		result.membranes.add( membrane );

		createNonMembraneMasks( result, t );
	}

	private RandomAccessibleInterval< BitType > computeCellEdge(
			TranslocationResult result,
			RandomAccessibleInterval< R > intensity )
	{
		final RandomAccessibleInterval< R > smoothed = Algorithms.median(
				intensity, blurRadius, opService );

		result.intensities.add( smoothed );

		final RandomAccessibleInterval< R > erode = Algorithms.erode(
				smoothed, erosionRadius );

		final RandomAccessibleInterval< R > gradient = Utils.createEmptyCopy( erode );

		LoopBuilder.setImages( gradient, smoothed, erode ).forEachPixel( ( g, i, e ) ->
				{
					g.setReal( i.getRealDouble() - e.getRealDouble() );
				}
		);

		result.gradients.add( gradient );

		RandomAccessibleInterval< BitType > binaryGradient =
				ArrayImgs2.img( intensity, new BitType(  ) );

		opService.threshold().isoData(
				Views.iterable( binaryGradient ),
				Views.iterable( gradient ) );

		Regions.onlyKeepLargestRegion( binaryGradient, FOUR_CONNECTED );

		// morphological smoothing for helping the thinning
		final RandomAccessibleInterval< BitType > openedBinaryGradient
				= Algorithms.dilate( Algorithms.erode( binaryGradient, 1 ), 1 );

		result.binaryGradients.add( openedBinaryGradient );

		return Skeletons.skeleton( openedBinaryGradient, opService );
	}

	private void measureMembraneIntensity( TranslocationResult result, int t )
	{
		final Cursor< BitType > cursor = Views.iterable(
				(RandomAccessibleInterval) result.membranes.get( t ) ).cursor();
		final RandomAccess< R > access = inputs.get( t ).randomAccess();

		double sum = 0.0;
		int n = 0;
		while ( cursor.hasNext() )
		{
			if ( cursor.next().get() )
			{
				access.setPosition( cursor );
				sum += access.get().getRealDouble();
				n++;
			}
		}

		result.membraneIntensities.add( sum / n );
	}

	private void createNonMembraneMasks( TranslocationResult result, int t )
	{
		final RandomAccessibleInterval< BitType > membraneSkeleton =
				( RandomAccessibleInterval< BitType > ) result.cellEdges.get( t );

		final RandomAccessibleInterval< BitType > nonMembranePixels =
				Converters.convert(
						// dilate to exclude membrane blur
						Algorithms.dilate( membraneSkeleton, 2 * resolutionBlurWidthInPixel ),
						( i, o ) -> o.set( !i.get() ),
						new BitType() );

		final ArrayList< RegionAndSize > sizeSortedRegions =
				Regions.getSizeSortedRegions(
						nonMembranePixels,
						FOUR_CONNECTED );

		final ArrayList< RandomAccessibleInterval< BitType > > nonMembraneMasks = new ArrayList<>();
		for ( int region = 0; region < Math.min( 2, sizeSortedRegions.size() ); region++ )
		{
			final RandomAccessibleInterval< BitType > nonMembraneMask = Utils.createEmptyCopy( membraneSkeleton );
			Regions.drawRegionInMask( sizeSortedRegions.get( region ).getRegion(), nonMembraneMask );
			nonMembraneMasks.add( nonMembraneMask );
		}

		result.nonMembraneMasks.add( nonMembraneMasks );
	}

	private void measureIntensities(
			TranslocationResult< ? > result,
			int t )
	{
		measureMembraneIntensity( result, t );

		measureNonMembraneIntensities( result, t );
	}

	private void measureNonMembraneIntensities(
			TranslocationResult< ? > result,
			int t )
	{
		final RandomAccessibleInterval< R > intensity = inputs.get( t );

		final ArrayList< RandomAccessibleInterval< BitType > > nonMembraneMasks =
				result.nonMembraneMasks.get( t );

		final ArrayList< Double > intensities = new ArrayList<>();

		if ( nonMembraneMasks.size() < 2 )
		{
			Logger.log( "Error: region " + ( results.size() + 1 ) + ", frame " + ( t + 1) );
			intensities.add( SEGMENTATION_ERROR );
			intensities.add( SEGMENTATION_ERROR );
		}
		else
		{
			for ( int region = 0; region < 2; region++ )
			{
				Measurements.sumIntensity( intensity, nonMembraneMasks.get( region )  );
				intensities.add( Measurements.sumIntensity( intensity, nonMembraneMasks.get( region )  ) );
			}

			Collections.sort( intensities );
		}

		result.dimmerIntensities.add( intensities.get( 0 ) );
		result.brighterIntensities.add( intensities.get( 1 ) );
	}

	private double getMean( ArrayList< RegionAndSize > regionsAndSizes,
							RandomAccess< R > intensityAccess,
							int i )
	{
		final LabelRegionCursor cursor = regionsAndSizes.get( i ).getRegion().cursor();

		double sum = 0.0;
		int n = 0;

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			intensityAccess.setPosition( cursor );
			sum += intensityAccess.get().getRealDouble();
			n++;
		}

		return sum / n;
	}

	private void computeIntensities( RandomAccessibleInterval< R > image, TranslocationResult result )
	{
		final RandomAccessibleInterval< BitType > mask = createCellMask( image );
		result.cellMasks.add( mask );

		final RandomAccessibleInterval< IntType > distances = computeDistanceMap( mask );

		result.dimmerIntensities.add( computeMeanIntensityWithinDistanceRange(
				image,
				distances,
				0,
				0 ) );

		final CoordinateAndValue membraneDistanceAndIntensity =
				computeMembraneIntensity( image, distances );

		result.membraneIntensities.add( membraneDistanceAndIntensity.value );

		result.brighterIntensities.add( computeMeanIntensityWithinDistanceRange(
				image,
				distances,
				(int) membraneDistanceAndIntensity.coordinate + resolutionBlurWidthInPixel,
				Integer.MAX_VALUE ) );
	}

	private CoordinateAndValue computeMembraneIntensity(
			RandomAccessibleInterval< R > image,
			RandomAccessibleInterval< IntType > distances )
	{
		double maximalGradient = 0;
		double membraneIntensity = 0;
		double distanceAtMaximalGradient = 0;

		ArrayList< Double > intensities = new ArrayList<>(  );

		final double dMax = Algorithms.getMaximumValue( distances );
		for ( int d = 0; d <= dMax; d++ )
		{
			double intensity = computeMeanIntensityWithinDistanceRange( image, distances, d, d );

			intensities.add( intensity );

			int previousDistance = d > resolutionBlurWidthInPixel ? d - resolutionBlurWidthInPixel : 0;

			double gradient = intensity - intensities.get( previousDistance );

			if ( gradient > maximalGradient )
			{
				maximalGradient = gradient;
				distanceAtMaximalGradient = d;
				membraneIntensity = intensity;
			}

		}

		return new CoordinateAndValue( distanceAtMaximalGradient, membraneIntensity );
	}

	private RandomAccessibleInterval< IntType > computeDistanceMap( RandomAccessibleInterval< BitType > mask )
	{
		final RandomAccessibleInterval< DoubleType > squaredDistances = Algorithms.computeSquaredDistances( mask );

		return Converters.convert(
				squaredDistances,
				( i, o ) -> o.set( (int) Math.sqrt( i.get() ) ),
				new IntType() );
	}

	private double computeMeanIntensityWithinDistanceRange(
			final RandomAccessibleInterval< R > image,
			final RandomAccessibleInterval< IntType > distances,
			final int dMin,
			final int dMax )
	{
		final Cursor< IntType > cursor = Views.iterable( distances ).cursor();
		final RandomAccess< R > access = image.randomAccess();

		double sum = 0.0;
		int n = 0;

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			if ( cursor.get().getInteger() >= dMin && cursor.get().getInteger() <= dMax )
			{
				access.setPosition( cursor );
				sum += access.get().getRealDouble();
				n++;
			}
		}

		if ( n == 0 )
		{
			final double maximumValue = Algorithms.getMaximumValue( distances );
			int a = 1;
		}

		return sum / n;

	}

	private RandomAccessibleInterval< BitType > createCellMask( RandomAccessibleInterval< R > image )
	{
		final SignalOverBackgroundSegmenter< R > segmenter =
				new SignalOverBackgroundSegmenter<>(
						image,
						signalToNoise,
						minimalRegionSize );

		RandomAccessibleInterval< BitType > mask = segmenter.createMask();

		mask = Algorithms.close( mask, closingRadius );

		return mask;
	}

}
