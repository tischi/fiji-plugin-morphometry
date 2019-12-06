package explore;

import de.embl.cba.morphometry.skeleton.Skeletons;
import net.imagej.ImageJ;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;

public class SkeletonsLongestBranchExample
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		RandomAccessibleInterval< BitType > skeleton = ArrayImgs.bits( 100, 100 );

		drawSkeleton02( skeleton );

		ImageJFunctions.show( skeleton, "skeleton" );

		final RandomAccessibleInterval< BitType > branchpoints = Skeletons.branchPoints( skeleton );

		ImageJFunctions.show( branchpoints, "branchPoints" );

		final RandomAccessibleInterval< BitType > longestBranch = Skeletons.longestBranch( skeleton );

		ImageJFunctions.show( longestBranch, "skeleton" );

	}

	public static void drawSkeleton01( RandomAccessibleInterval< BitType > skeleton )
	{
		final RandomAccess< BitType > access = skeleton.randomAccess();

		access.setPosition( 30, 1 );
		for ( int x = 10; x < 90; x++ )
		{
			access.setPosition( x, 0 );
			access.get().set( true );
		}

		access.setPosition( 20, 0 );
		for ( int y = 20; y < 40; y++ )
		{
			access.setPosition( y, 1 );
			access.get().set( true );
		}
	}

	public static void drawSkeleton02( RandomAccessibleInterval< BitType > skeleton )
	{
		final RandomAccess< BitType > access = skeleton.randomAccess();

		access.setPosition( 30, 1 );
		for ( int x = 10; x < 50; x++ )
		{
			access.setPosition( x, 0 );
			access.get().set( true );
		}

		access.setPosition( 31, 1 );
		for ( int x = 50; x < 70; x++ )
		{
			access.setPosition( x, 0 );
			access.get().set( true );
		}


	}
}
