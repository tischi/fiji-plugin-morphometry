package explore;

import java.awt.*;

public class GetColorFromName
{
	public static void main( String[] args ) throws NoSuchFieldException, IllegalAccessException
	{
		final Color white = (Color)Color.class.getField("WHITE").get(null);

		final Color green = Color.getColor( "GREEN" );
		final Color magenta = Color.getColor( "Magenta" );
	}
}
