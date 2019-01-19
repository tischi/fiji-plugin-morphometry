package de.embl.cba.morphometry.skeleton;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.regions.Regions;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.morphology.table2d.Branchpoints;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.logic.BitType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public class Skeletons
{
	public static RandomAccessibleInterval< BitType > longestBranch(
			RandomAccessibleInterval< BitType > skeleton
	)
	{

		final RandomAccessibleInterval< BitType > longestBranch = Utils.copyAsArrayImg( skeleton );

		removeBranchpoints( longestBranch );

		Regions.onlyKeepLargestRegion( longestBranch,
				ConnectedComponents.StructuringElement.EIGHT_CONNECTED );

		return longestBranch;
	}

	public static RandomAccessibleInterval< BitType > branchpoints(
			RandomAccessibleInterval< BitType > skeleton )
	{
		RandomAccessibleInterval< BitType > branchpoints =
				ArrayImgs.bits( Intervals.dimensionsAsLongArray( skeleton ) );

		Branchpoints.branchpoints(
				Views.extendBorder( Views.zeroMin( skeleton ) ),
				Views.flatIterable( branchpoints ) );

		Views.translate( branchpoints, Intervals.minAsLongArray( skeleton ) );

		return branchpoints;
	}

	static void removeBranchpoints( RandomAccessibleInterval< BitType > skeleton )
	{
		final RandomAccessibleInterval< BitType > branchpoints = branchpoints( skeleton );

		LoopBuilder.setImages( branchpoints, skeleton ).forEachPixel( ( b, s ) ->
		{
			if ( b.get() ) s.set( false );
		});
	}


	public static RandomAccessibleInterval< BitType > skeleton(
			RandomAccessibleInterval< BitType > input,
			OpService opService )
	{
		final RandomAccessibleInterval< BitType > thin = Utils.createEmptyCopy( input );
		opService.morphology().thinGuoHall( thin, input );
		return thin;
	}
}
