package de.embl.cba.morphometry.splitting;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.SyncWindowsHack;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaTrackingSettings;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.frame.SyncWindows;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;

public class TrackingSplitter< T extends RealType< T > & NativeType< T > >
{

	final ArrayList< RandomAccessibleInterval< BitType > > masks;
	final ArrayList< RandomAccessibleInterval< T > > intensities;

	private int nextId;

	private ArrayList< RandomAccessibleInterval< BitType > > splitMasks;
	final MicrogliaTrackingSettings settings;

	public TrackingSplitter( ArrayList< RandomAccessibleInterval< BitType > > masks,
							 ArrayList< RandomAccessibleInterval< T > > intensities,
							 MicrogliaTrackingSettings settings )
	{
		this.masks = masks;
		this.intensities = intensities;
		this.settings = settings;
	}


	public void run()
	{

		int tMin = 0;  // at this point the movie is already cropped in time, such that we can process the full movie
		int tMax = masks.size() - 1;
		int t = tMin;

		splitMasks = new ArrayList<>( );


		/**
		 * Process first time-point
		 */

		Utils.log( "\nRunning ShapeAndIntensitySplitter on frame " + settings.tMin );
		final ShapeAndIntensitySplitter splitter = new ShapeAndIntensitySplitter( masks.get( t ), intensities.get( t ), settings );
		splitter.run();

		/**
		 * Allow user to manually correct
		 */

		if ( settings.manualSegmentationCorrectionOfFirstFrame )
		{
			IJ.run("Brightness/Contrast...");
			splitMasks.add( getManuallyCorrectedMask( splitter.getSplitMask(), t ) );
		}
		else
		{
			splitMasks.add( splitter.getSplitMask() );
		}

		nextId = Utils.getNumObjects( splitMasks.get( tMin ) );

		boolean showSplittingAttempts = false;


		/**
		 * Process subsequent time-points
		 */

		RandomAccessibleInterval< IntType > previousLabeling = Utils.asImgLabeling( splitMasks.get( tMin ) ).getSource();

		for ( t = tMin + 1; t <= tMax; ++t )
		{

			Utils.log( "\nProcessing frame " + ( t + 1 ) );

			final ImgLabeling< Integer, IntType > currentImgLabeling = Utils.asImgLabeling( masks.get( t ) );
			RandomAccessibleInterval< IntType > currentLabeling = currentImgLabeling.getSource();

			HashMap< Integer, ArrayList< Integer > > overlappingObjectsLabelsMap = getOverlappingObjectLabelsMap( t, previousLabeling, currentImgLabeling, currentLabeling );

			RandomAccessibleInterval< BitType > splitMask = Utils.copyAsArrayImg( masks.get( t ) );

			Algorithms.splitCurrentObjectsBasedOnOverlapWithPreviousObjects(
					splitMask,
					overlappingObjectsLabelsMap,
					currentImgLabeling,
					intensities.get( t ),
					previousLabeling,
					( long ) ( settings.minimalTrackingSplittingObjectArea / Math.pow( settings.workingVoxelSize, splitMask.numDimensions() ) ),
					( int ) ( settings.minimalObjectCenterDistance / settings.workingVoxelSize ),
					settings.opService,
					false);

			if ( settings.manualSegmentationCorrectionOfAllFrames )
			{
				splitMask = getManuallyCorrectedMask( splitMask, t );
			}

			splitMasks.add( splitMask );

			previousLabeling = Utils.asImgLabeling( splitMask ).getSource();

		}
	}

	public RandomAccessibleInterval< BitType > getManuallyCorrectedMask( RandomAccessibleInterval< BitType > mask, int t )
	{
		final ImagePlus intensitiesImp = Utils.createIJ1Movie( intensities, "intensities" );
		intensitiesImp.show();
		intensitiesImp.setT( t + 1 );
		intensitiesImp.updateImage();
		//intensitiesImp.updateAndDraw();

		final RandomAccessibleInterval< IntType > indexImg = Utils.asImgLabeling( mask ).getIndexImg();
		ImagePlus labelImagePlus = Utils.asLabelImagePlus( indexImg );
		labelImagePlus.setTitle( "Label mask of frame " + ( t + 1 ) );
		labelImagePlus.show();
		IJ.run( labelImagePlus, "Enhance Contrast", "saturated=0.35");

		IJ.wait( 100 );
		final SyncWindowsHack syncWindows = new SyncWindowsHack();
		syncWindows.syncAll();

		final NonBlockingGenericDialog gd = new NonBlockingGenericDialog( "Manual correction" );
		gd.addMessage( "Please correct segmentation of frame " + ( t + 1 )
				+ " and click [Ok] when you are done.\n"
				+ "Please click [Cancel] in order to skip manual correction of all proceeding frames." );
		gd.showDialog();


		if ( gd.wasCanceled() )
		{
			settings.manualSegmentationCorrectionOfAllFrames = false;
		}

		syncWindows.close();

		intensitiesImp.close();
		labelImagePlus.hide();

		return Utils.asMask( (RandomAccessibleInterval) ImageJFunctions.wrapReal( labelImagePlus ) );
	}

	public HashMap< Integer, ArrayList< Integer > > getOverlappingObjectLabelsMap( int t, RandomAccessibleInterval< IntType > previousLabeling, ImgLabeling< Integer, IntType > currentImgLabeling, RandomAccessibleInterval< IntType > currentLabeling )
	{
		HashMap< Integer, ArrayList< Integer > > overlappingObjectsLabelsMap = new HashMap<>(  );

		LabelRegions< Integer > labelRegions = new LabelRegions( currentImgLabeling );

		for ( LabelRegion< Integer > region : labelRegions )
		{
			final HashMap< Integer, Long > overlapSizes = computeOverlaps( previousLabeling, region );

			if ( overlapSizes.size() == 0 )
			{
				overlappingObjectsLabelsMap.put( region.getLabel(), new ArrayList<>() );
			}
			else
			{
				if ( overlapSizes.size() > 1 )
				{
					Utils.log( "Object at "
							+ " x = " + ( int ) region.getCenterOfMass().getDoublePosition( 0 )
							+ " y = " + ( int ) region.getCenterOfMass().getDoublePosition( 1 )
							+ " overlaps with " + overlapSizes.size() + " objects in previous frame." );
				}

				final ArrayList< Integer > trulyOverlappingObjectLabels
						= getTrulyOverlappingObjectLabels(
								intensities.get( t - 1 ),
								intensities.get( t ),
								previousLabeling,
								currentLabeling,
								region,
								overlapSizes );

				overlappingObjectsLabelsMap.put( region.getLabel(), trulyOverlappingObjectLabels );

			}

		}
		return overlappingObjectsLabelsMap;
	}


	public ArrayList< Integer > getTrulyOverlappingObjectLabels( RandomAccessibleInterval< T > previousIntensityImage,
												RandomAccessibleInterval< T > currentIntensityImage,
												RandomAccessibleInterval< IntType > previousLabeling,
												RandomAccessibleInterval< IntType > currentLabeling,
												LabelRegion< Integer > currentObjectRegion,
												HashMap< Integer, Long > overlapSizes )
	{

		final ArrayList< Integer > trulyOverlappingObjectLabels = new ArrayList<>();

		if ( overlapSizes.keySet().size() == 1 )
		{
			trulyOverlappingObjectLabels.add( overlapSizes.keySet().iterator().next() );
		}
		else
		{
			for ( int previousLabel : overlapSizes.keySet() )
			{
				final long previousObjectSize = Measurements.measureSize( previousLabeling, previousLabel );

				final double overlapFraction = 1.0 * overlapSizes.get( previousLabel ).longValue() / previousObjectSize;

				// Utils.log( "Previous object size / overlap size: " + overlapFraction );

				if ( overlapFraction >= settings.minimalOverlapFraction )
				{
					trulyOverlappingObjectLabels.add( previousLabel );
				}
			}

			// Utils.log( "Corrected number of overlapping objects: " + trulyOverlappingObjectLabels.size() );

		}

		return trulyOverlappingObjectLabels;
	}

	public boolean isReallyTwoObjects( RandomAccessibleInterval< T > previousIntensityImage,
									   RandomAccessibleInterval< T > currentIntensityImage,
									   RandomAccessibleInterval< IntType > previousLabeling,
									   RandomAccessibleInterval< IntType > currentLabeling,
									   LabelRegion< Integer > currentObjectRegion,
									   HashMap< Integer, Long > previousSizes )
	{

		boolean splitObjects =  true;

		final double currentObjectIntensity = Measurements.measureBgCorrectedSumIntensity(
				currentLabeling,
				currentObjectRegion.getLabel(),
				currentIntensityImage );

		Utils.log( "Object intensity: " + (long) currentObjectIntensity );

		final HashMap< Integer, Double > previousIntensities = new HashMap<>();

		for ( int previousLabel : previousSizes.keySet() )
		{
			final double previousObjectIntensity = Measurements.measureBgCorrectedSumIntensity(
					previousLabeling,
					previousLabel,
					previousIntensityImage );

			previousIntensities.put( previousLabel , previousObjectIntensity );

			final double overlapFraction = 1.0 * previousSizes.get( previousLabel ).longValue() / currentObjectRegion.size();
			Utils.log( "Previous object intensity: " + (long) previousObjectIntensity );
			Utils.log( "Overlap pixel fraction: " + overlapFraction );

			if ( overlapFraction < settings.minimalOverlapFraction ) splitObjects = false;

		}

		double previousIntensitySum = getPreviousIntensitySum( previousIntensities );

		final double sumIntensityRatio = currentObjectIntensity / previousIntensitySum;

		Utils.log( "Intensity ratio: " + sumIntensityRatio );

		if ( sumIntensityRatio < settings.minimalSumIntensityRatio ) splitObjects = false;
		if ( sumIntensityRatio > settings.maximalSumIntensityRatio ) splitObjects = false;

		Utils.log( "Split objects: " + splitObjects );

		return splitObjects;
	}

	private double getPreviousIntensitySum( HashMap< Integer, Double > previousIntensities )
	{
		double previousIntensitySum = 0;

		for ( double intensity : previousIntensities.values() )
		{
			previousIntensitySum += intensity;

		}
		return previousIntensitySum;
	}

	public ArrayList< RandomAccessibleInterval< BitType > > getSplitMasks()
	{
		return splitMasks;
	}



	public ArrayList< RandomAccessibleInterval< IntType > > initUpdatedLabelings( int tMin )
	{
//		final ArrayList< RandomAccessibleInterval< IntType > > updatedLabelings = new ArrayList<>();
//		updatedLabelings.add( masks.get( tMin ).getSource() );
//		return updatedLabelings;
		return null;
	}

	public int computeObjectId( HashMap< Integer, Long > overlaps )
	{
		int objectId;

		if( overlaps.size() == 0 )
		{
			objectId = nextId++;
		}
		else
		{
			objectId = getMaxOverlapLabel( overlaps );
		}

		return objectId;
	}

	public int getMaxOverlapLabel( HashMap< Integer, Long > overlaps )
	{
		int maxOverlapLabel = 0;
		long maxOverlapCount = Long.MIN_VALUE;

		for ( Integer label : overlaps.keySet() )
		{
			final long overlapCount = overlaps.get( label ).longValue();

			if ( overlapCount > maxOverlapCount )
			{
				maxOverlapCount = overlapCount;
				maxOverlapLabel = label;
			}
		}
		return maxOverlapLabel;
	}

	public HashMap< Integer, Long > computeOverlaps(
			RandomAccessibleInterval< IntType > previousLabeling,
			LabelRegion region )
	{
		final HashMap< Integer, Long > overlaps = new HashMap<>();

		final RandomAccess< IntType > previousLabelsAccess = previousLabeling.randomAccess();
		final LabelRegionCursor cursor = region.cursor();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			previousLabelsAccess.setPosition( cursor );
			final int previousLabel = previousLabelsAccess.get().getInteger();

			if ( previousLabel != 0 )
			{
				addOverlap( overlaps, previousLabel );
			}
		}

		return overlaps;
	}


	public void addOverlap( HashMap< Integer, Long > overlaps, int integer )
	{
		if ( overlaps.keySet().contains( integer ) )
		{
			overlaps.put( integer,  overlaps.get( integer ).longValue() + 1 );
		}
		else
		{
			overlaps.put( integer,  1L );
		}
	}


}
