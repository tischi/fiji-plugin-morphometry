package de.embl.cba.morphometry.translocation;

import de.embl.cba.tables.models.ColumnClassAwareTableModel;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;

import javax.swing.*;
import java.util.ArrayList;

public class TranslocationResult < T extends RealType< T > & NativeType< T > >
{
	public static final String REGION_INDEX = "RegionIndex";
	public static final String TIMEPOINT = "Timepoint";
	public static final String TRANSLOCATION = "Translocation";
	public static final String REGION_CENTER_X = "RegionCenterX";
	public static final String REGION_CENTER_Y = "RegionCenterY";
	final public ArrayList< RandomAccessibleInterval< T > > cellMasks;
	final public ArrayList< RandomAccessibleInterval< T > > gradients;
	final public ArrayList< RandomAccessibleInterval< BitType > > membraneMasks;
	final public ArrayList< RandomAccessibleInterval< BitType > > insideOutsideMasks;

	final public ArrayList< Double > outsideIntensities;
	final public ArrayList< Double > membraneIntensities;
	final public ArrayList< Double > insideIntensities;
	final public ArrayList< Double > translocations;
	public double regionCenterX;
	public double regionCenterY;

	public TranslocationResult( )
	{
		this.cellMasks = new ArrayList<>(  );
		outsideIntensities = new ArrayList<Double>();
		membraneIntensities = new ArrayList<Double>();
		insideIntensities = new ArrayList<Double>();
		translocations = new ArrayList<Double>();
		gradients = new ArrayList<RandomAccessibleInterval<T>>();
		membraneMasks = new ArrayList<RandomAccessibleInterval<BitType>>();
		insideOutsideMasks = new ArrayList<RandomAccessibleInterval<BitType>>();
	}

	public static JTable resultsAsTable( ArrayList< TranslocationResult > results )
	{
		final int numRegions = results.size();
		int numTimepoints = results.get( 0 ).translocations.size();

		String[] columnNames = new String[]{
				REGION_INDEX,
				REGION_CENTER_X,
				REGION_CENTER_Y,
				TIMEPOINT,
				TRANSLOCATION
		};

		final ColumnClassAwareTableModel model =
				new ColumnClassAwareTableModel(
						numRegions * numTimepoints,
						columnNames.length
				);

		model.setColumnIdentifiers( columnNames );

		model.addRow( new Object[ numRegions * numTimepoints ] );

		for ( int r = 0; r < numRegions; r++ )
		{

			for ( int t = 0; t < numTimepoints; t++ )
			{

				final int row = t + numTimepoints * r;

				model.setValueAt(
						r + 1,
						row,
						model.findColumn( REGION_INDEX ));

				model.setValueAt(
						results.get( r ).regionCenterX,
						row,
						model.findColumn( REGION_CENTER_X ));

				model.setValueAt(
						results.get( r ).regionCenterY,
						row,
						model.findColumn( REGION_CENTER_Y ));

				model.setValueAt(
						t,
						row,
						model.findColumn( TIMEPOINT ));

				model.setValueAt(
						results.get( r ).translocations.get( t ),
						row,
						model.findColumn( TRANSLOCATION ) );

			}
		}

		model.refreshColumnClasses();

		return new JTable( model );
	}
}
