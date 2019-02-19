package de.embl.cba.morphometry.geometry.ellipsoids;

public class EllipsoidMLJ
{
	public static final int PHI = 0, THETA = 1, PSI = 2;

	public double[] center = new double[ 3 ];
	public double[] radii = new double[ 3 ];
	public double[] eulerAnglesInDegrees = new double[ 3 ];

	@Override
	public String toString()
	{
		String s = "";
		s += "\nMLJ ellipsoid angles (Phi, Theta, Psi):";
		s += "\nPhi: " + eulerAnglesInDegrees[ PHI ];
		s += "\nTheta: " + eulerAnglesInDegrees[ THETA ];
		s += "\nPsi: " + eulerAnglesInDegrees[ PSI ];
		return s;
	}
}
