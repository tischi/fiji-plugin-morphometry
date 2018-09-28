package de.embl.cba.morphometry.splitting;

import de.embl.cba.morphometry.measurements.ObjectMeasurements;
import de.embl.cba.morphometry.microglia.MicrogliaSettings;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.HashMap;
import java.util.Map;

public class SplittingUtils
{

	public static HashMap< Integer, Integer > getNumObjectsFromSkeleton( ImgLabeling< Integer, IntType > imgLabeling, RandomAccessibleInterval< BitType > skeleton, MicrogliaSettings settings )
	{
		HashMap< Integer, Map< String, Object > > skeletonMeasurements = new HashMap<>();
		ObjectMeasurements.measureSumIntensities( skeletonMeasurements, imgLabeling, skeleton, "skeleton" );
		HashMap< Integer, Integer > numObjects = new HashMap<>();
		for ( int label : skeletonMeasurements.keySet() )
		{
			final double skeletonLength = settings.workingVoxelSize * (long) skeletonMeasurements.get( label ).get( ObjectMeasurements.SUM_INTENSITY + "_skeleton" );
			int n = (int) ( Math.ceil( skeletonLength / settings.skeletonMaxLength ) );
			numObjects.put( label, n );
		}
		return numObjects;
	}
}
