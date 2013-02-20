/* AUSLAGERUNG DER TESTROUTINEN
 * Funktionenheader nach wie vor in osmosisSetup.h
 * um Zugriff auf die globalen Klassenvariablen zu haben
 */

namespace DROPS {

/// Compute Mean Drop Concentration, i.e. mean
/// concentration in second phase (negative sign):
/// integral over concentration / volume of second phase
/*double OsmosisP1CL::MeanDropConcentration(bool absolute)
{
    VecDescCL cn (&idx);
    GetSolutionOnPart(cn, false, false);
    double absdet;
    InterfaceTetraCL patch;

    double c_avrg= 0., Volume= 0.;
    LocalP2CL<double> ones( 1.);
    const Uint lvl= ct.GetLevel();
    DROPS_FOR_TRIANG_TETRA( MG_, lvl, it) {
        LocalP1CL<> lp1_cn( *it, cn, Bnd_);
        LocalP2CL<> lp2_cn( lp1_cn );
        absdet= std::abs( it->GetVolume()*6.);
        patch.Init( *it, lset_,0.);
        for (int ch=0; ch<8; ++ch)
        {
            // compute volume and concentration
            patch.ComputeCutForChild(ch);
            Volume+= patch.quad( ones, absdet, false);
            c_avrg+= patch.quad( lp2_cn, absdet, false);
        }
    }
    if (!absolute)
    	c_avrg/= Volume;
    return c_avrg;
}*/


/*double  AnalyticSolutionCase1(const DROPS::Point3DCL& p, double t)
{
	static double r0 = 1.0;
	static double dn = P.get<double>("Osmosis.Diffusivity");
	static double vn = P.get<double>("Osmosis.GrowVelocity");
	static DROPS::Point3DCL drop = P.get<DROPS::Point3DCL>("Exp.PosDrop");
	double r2 = (p[0]-drop[0])*(p[0]-drop[0])
		+ (p[1]-drop[1])*(p[1]-drop[1]);
	//		  + (p[2]-drop[2])*(p[2]-drop[2]);
	//EXPANSION / CONTRACTION
	//return vn/dn * r2 - 2*(r0+vn*t) - vn/dn * (r0+vn*t)*(r0+vn*t);
	//STATIC
	return cos(r2 - 1) + t*t;
}*/

/*double  AnalyticVolumeCase1(double& t)
{
	return 4./3.*M_PI*(1.+t)*(1.+t)*(1.+t);
}*/

/*double  AnalyticVolumeCase1(double fac, double& t)
{
	return fac*4./3.*M_PI*(2.-t)*(2.-t)*(2.-t);
}*/

void OsmosisP1CL::SolutionErrorCase1(double& con, double& vol, double )
{
    double absdet;
    InterfaceTetraCL patch;

    con = 0.0;
    vol = 0.0;
    LocalP2CL<double> ones( 1.);
    const Uint lvl= conc.GetLevel();
    DROPS_FOR_TRIANG_TETRA( MG_, lvl, it) {
        LocalP1CL<> lp1_cn( *it, conc, Bnd_);
        LocalP2CL<> lp2_cn( lp1_cn );
        absdet= std::abs( it->GetVolume()*6.);
        patch.Init( *it, lset_.Phi,0.);
        for (int ch=0; ch<8; ++ch)
        {
            // compute volume and concentration
            patch.ComputeCutForChild(ch);
            //Volume += patch.quad( ones, absdet, false);
            //c_avrg += std::abs(patch.quad( lp2_cn, absdet, false));
            con  += std::abs(patch.quad( lp2_cn, absdet, false));
			vol  += std::abs(patch.quad( ones, absdet, false));
        }
    }
}


Point3DCL TestFunctionTemp(const Point3DCL& p, double )
{
	static Point3DCL pos = P.get<Point3DCL>("Exp.PosDrop");
	return 2*(p-pos);
}

double TestFunctionConc(const Point3DCL& p, double )
{
	//static Point3DCL pos = P.get<Point3DCL>("Exp.PosDrop");
	return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
}

Point3DCL AnalyticNormal(const Point3DCL& p, double )
{
	static Point3DCL pos = P.get<Point3DCL>("Exp.PosDrop");
	return 2*(p-pos);
}

void OsmosisP1CL::SolutionErrorCase2(double& verr, double& aerr, double )
{
	std::cout << "Calculating Error...\n";
	double det;
	InterfaceTriangleCL tri;

	double analyticArea = 4*M_PI*P.get<Point3DCL>("Exp.RadDrop")[0];

	verr = 0.;
	aerr = 0.;
	double v = 0.;
	LocalP2CL<double> ones( 1.);

	DROPS_FOR_TRIANG_TETRA( MG_, MG_.GetLastLevel(), sit)
	{

		//compute triangles that form the approximation of the interface
		tri.Init(*sit, lset_.Phi, 0.);

		if (!tri.Intersects())
			continue;

		//Boundary probably not right!
		LocalP1CL<Point3DCL>	solP1(*sit, Vn_, Bnd_v_);
		LocalP2CL<Point3DCL> 	solP2(solP1);

		//loop over the 8 sub-tetras (childs) that were constructed
		for (int ch= 0; ch < 8; ++ch)
		{
			//continue if there exists no intersection with the interface on the child
			if (!tri.ComputeForChild( ch))
				continue;
			//switched to child now

			det = tri.GetAbsDet();

			//loop over triangles that form the patch (maximum of 2)
			for (int t= 0; t < tri.GetNumTriangles(); ++t)
			{

				if (t == 1) //Quadrilateral
					det *= tri.GetAreaFrac();
				if (t > 1)
					std::cout << "Quadrilateral divided in > 2 triangles, should not happen\n";

				Quad5_2DCL<Point3DCL> 	analytic(*sit, &tri.GetBary(t), &TestFunctionTemp);
				Quad5_2DCL<double> 		area(ones, &tri.GetBary(t));
				Quad5_2DCL<Point3DCL> 	sol(solP2, &tri.GetBary(t));
				Quad5_2DCL<Point3DCL>	qerr2((analytic-sol)*(analytic-sol));
				//first we need quadrature objects for the ansatz functions

				verr += qerr2.quad(det).norm();
				aerr += area.quad(det);
				v += analytic.quad(det).norm();

			}
		}
	}

	aerr -= analyticArea;
	verr = std::sqrt(verr)/v;
    std::cout << "Area error  :  " << std::setprecision(12) << aerr <<"\n";
    std::cout << "Normal error:  " << std::setprecision(12) << verr <<"\n";
}

void OsmosisP1CL::SolutionErrorCase3(double& verr, double& aerr, double )
{
	std::cout << "Calculating Error...\n";
	double det;
	InterfaceTriangleCL tri;

	verr = 0.;
	aerr = 0.;
	double v = 0.0;

	DROPS_FOR_TRIANG_TETRA( MG_, MG_.GetLastLevel(), sit)
	{

		//compute triangles that form the approximation of the interface
		tri.Init(*sit, lset_.Phi, 0.);

		if (!tri.Intersects())
			continue;

		//Boundary probably not right!
		LocalP1CL<Point3DCL>	solP1(*sit, Vn_, Bnd_v_);
		LocalP2CL<Point3DCL> 	solP2(solP1);

		//loop over the 8 sub-tetras (childs) that were constructed
		for (int ch= 0; ch < 8; ++ch)
		{
			//continue if there exists no intersection with the interface on the child
			if (!tri.ComputeForChild( ch))
				continue;
			//switched to child now

			det = tri.GetAbsDet();

			//loop over triangles that form the patch (maximum of 2)
			for (int t= 0; t < tri.GetNumTriangles(); ++t)
			{

				if (t == 1) //Quadrilateral
					det *= tri.GetAreaFrac();
				if (t > 1)
					std::cout << "Quadrilateral divided in > 2 triangles, should not happen\n";

				Quad5_2DCL<Point3DCL>	normal(*sit, &tri.GetBary(t), &AnalyticNormal);
				Quad5_2DCL<Point3DCL> 	sol(solP2, &tri.GetBary(t));
				Quad5_2DCL<Point3DCL>	qerr2((normal-sol)*(normal-sol));
				//first we need quadrature objects for the ansatz functions

				verr += qerr2.quad(det).norm();
				v += normal.quad(det).norm();

			}
		}
	}
	verr = std::sqrt(verr)/v;

    std::cout << "relative error:  " << std::setprecision(12) << verr <<"\n";
}

void OsmosisP1CL::SolutionErrorCase4(double& verr, double& aerr, double t)
{
	std::cout << "Calculating Error...\n";
	double det, absdet;
	InterfaceTriangleCL tri;
	InterfaceTetraCL patch;

	double analyticVolume = 4./3.*M_PI*pow(std::sqrt(4.-4.*t), 3);

	verr = 0.;
	aerr = 0.;
	double Volume = 0.;
	double Area = 0.;
	LocalP2CL<double> ones( 1.);

	DROPS_FOR_TRIANG_TETRA( MG_, MG_.GetLastLevel(), sit)
	{

		//compute triangles that form the approximation of the interface
		tri.Init(*sit, lset_.Phi, 0.);

		absdet= std::abs( sit->GetVolume()*6.);
		patch.Init( *sit, lset_.Phi,0.);

		//loop over the 8 sub-tetras (childs) that were constructed
		for (int ch= 0; ch < 8; ++ch)
		{
			patch.ComputeCutForChild(ch);
			Volume += patch.quad( ones, absdet, false);
		}
	}

	DROPS_FOR_TRIANG_TETRA( MG_, MG_.GetLastLevel(), sit)
	{

		//compute triangles that form the approximation of the interface
		tri.Init(*sit, lset_.Phi, 0.);

		//loop over the 8 sub-tetras (childs) that were constructed
		for (int ch= 0; ch < 8; ++ch)
		{
			//continue if there exists no intersection with the interface on the child
			if (!tri.ComputeForChild( ch))
				continue;
			//switched to child now


			det = tri.GetAbsDet();

			//loop over triangles that form the patch (maximum of 2)
			for (int t= 0; t < tri.GetNumTriangles(); ++t)
			{

				if (t == 1) //Quadrilateral
					det *= tri.GetAreaFrac();
				if (t > 1)
					std::cout << "Quadrilateral divided in > 2 triangles, should not happen\n";

				Quad5_2DCL<double> 		area(ones, &tri.GetBary(t));

				Area += area.quad(det);
			}
		}
	}


	aerr = pow(M_PI, 1./3.)*pow((6*Volume), 2./3.) / Area;
	verr = std::abs((Volume - analyticVolume)/analyticVolume);
    std::cout << "Sphericity:  " << std::setprecision(12) << aerr <<"\n";
    std::cout << "Area:  " << std::setprecision(12) << Area <<"\n";
    std::cout << "ana. Area:  " << std::setprecision(12) << 4.*M_PI*pow(std::sqrt(4.-4.*t), 2) <<"\n";
    std::cout << "Volume:  " << std::setprecision(12) << Volume <<"\n";
    std::cout << "ana. Volume:  " << std::setprecision(12) << analyticVolume <<"\n";
    std::cout << "Volume error:  " << std::setprecision(12) << verr <<"\n";
}

}
