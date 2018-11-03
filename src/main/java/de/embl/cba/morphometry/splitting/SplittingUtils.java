package de.embl.cba.morphometry.splitting;

import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.microglia.MicrogliaTrackingSettings;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.HashMap;
import java.util.Map;

public class SplittingUtils
{

	public static HashMap< Integer, Integer > getNumObjectsFromSkeleton( ImgLabeling< Integer, IntType > imgLabeling, RandomAccessibleInterval< BitType > skeleton, MicrogliaTrackingSettings settings )
	{
		HashMap< Integer, Map< String, Object > > skeletonMeasurements = new HashMap<>();
		Measurements.measureSumIntensities( skeletonMeasurements, imgLabeling, skeleton, "Skeleton" );
		HashMap< Integer, Integer > numObjects = new HashMap<>();
		for ( int label : skeletonMeasurements.keySet() )
		{
			final double skeletonLength = settings.workingVoxelSize * (long) skeletonMeasurements.get( label ).get( Measurements.SUM_INTENSITY + "_skeleton" );
			int n = (int) ( Math.ceil( skeletonLength / settings.skeletonMaxLength ) );
			numObjects.put( label, n );
		}
		return numObjects;
	}
}
