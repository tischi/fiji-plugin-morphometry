package de.embl.cba.morphometry.splitting;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaSegmentationAndTrackingSettings;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.HashMap;

import static de.embl.cba.morphometry.splitting.SplittingUtils.getNumObjectsFromSkeleton;


public class ShapeAndIntensitySplitter< T extends RealType< T > & NativeType< T > >
{

	final MicrogliaSegmentationAndTrackingSettings settings;
	final private RandomAccessibleInterval< BitType > mask;
	final private RandomAccessibleInterval< T > intensity;
	private RandomAccessibleInterval< BitType > splitMask;

	public ShapeAndIntensitySplitter( RandomAccessibleInterval< BitType > mask,
									  RandomAccessibleInterval< T > intensity,
									  MicrogliaSegmentationAndTrackingSettings settings )
	{
		this.mask = mask;
		this.intensity = intensity;
		this.settings = settings;
	}

	public void run()
	{

		/**
		 * Get objects
		 */

		final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( mask );

		/**
		 * Estimate number of objects from skeleton
		 */

		// TODO: implement skeleton per object such that one can do closing operations without joining neighboring objects

		RandomAccessibleInterval< BitType > skeleton = settings.opService.morphology().thinGuoHall( mask );

		HashMap< Integer, Integer > numObjectsPerRegion = getNumObjectsFromSkeleton( imgLabeling, skeleton, settings );

		splitMask = Utils.copyAsArrayImg( mask );

		Algorithms.splitTouchingObjects(
				imgLabeling,
				intensity,
				splitMask,
				numObjectsPerRegion,
				( int ) ( settings.minimalObjectCenterDistance / settings.workingVoxelSize ),
				( long ) ( settings.minimalObjectSize / Math.pow( settings.workingVoxelSize, intensity.numDimensions() ) ),
				( int ) ( settings.maximalWatershedLength / settings.workingVoxelSize ),
				settings.opService, false, false );

	}

	public RandomAccessibleInterval< BitType > getSplitMask()
	{
		return splitMask;
	}


}
