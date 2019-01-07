package de.embl.cba.morphometry.translocation;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class TranslocationResult < T extends RealType< T > & NativeType< T > >
{
	public ArrayList< RandomAccessibleInterval< T > > cellMasks;

	public TranslocationResult( )
	{
		this.cellMasks = cellMasks;
	}
}
