import net.imglib2.util.LinAlgHelpers;

import static java.lang.Math.acos;

public class TestRotationAngle
{
	public static void main( String[] args )
	{

		double[] xAxis = new double[]{ 1, 0, 0};

		double[] vector = new double[]{ -1, -1, 0};

		final double dot = LinAlgHelpers.dot( xAxis, vector );

		double angleInRadians = acos( dot / ( LinAlgHelpers.length( xAxis ) * LinAlgHelpers.length( vector ) ) );

		System.out.println( 180 / Math.PI *  angleInRadians );

	}
}
