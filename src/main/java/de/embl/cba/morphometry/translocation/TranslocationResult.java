package de.embl.cba.morphometry.translocation;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class TranslocationResult < T extends RealType< T > & NativeType< T > >
{
	final public ArrayList< RandomAccessibleInterval< T > > cellMasks;
	final public ArrayList< Double > outsideIntensities;
	final public ArrayList< Double > membraneIntensities;
	final public ArrayList< Double > insideIntensities;
	final public ArrayList< Double > translocation;

	public TranslocationResult( )
	{
		this.cellMasks = new ArrayList<>(  );
		outsideIntensities = new ArrayList<Double>();
		membraneIntensities = new ArrayList<Double>();
		insideIntensities = new ArrayList<Double>();
		translocation = new ArrayList<Double>();
	}
}
