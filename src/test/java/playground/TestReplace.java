package playground;

public class TestReplace
{
	public static void main( String[] args )
	{
		String a = "aaaaa-1.tif";

		final String replace = a.replace( ".tif", "" );

		System.out.println( replace );
	}
}
