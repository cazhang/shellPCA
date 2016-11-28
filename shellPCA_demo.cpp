#include <aol.h>
#include <parameterParser.h>
#include <eigenvectors.h>
#include <multiStreambuf.h>

#include "../../../finishedProjects/geodesicCalculusShells/averagingOps.h"

typedef double RType;
typedef om::TriMesh<RType> TriMeshType;
typedef aol::TriangMeshConfigurator<RType, TriMeshType, aol::CenterQuadrature <RType> > TriMeshCenterQuadConf;
typedef aol:: SparseMatrix<RType> SubMatrixType;

//! main
int main(int argc, char **argv)
{
	try{
		//! set path and paras
		//! to run, give two arguments
		if ( argc != 2 ) { 
			cerr << "Wrong number of input files." << endl;
			cerr << "Usage: ./shellPCA_demo <parameterFile>" << endl;
			return EXIT_FAILURE;
		}
		// read in file names of parameter files
		char parsername[1024];
		sprintf( parsername, "%s",  argv[1] );
		cerr << endl << "Reading parameters from '" << parsername << "'." << endl;
		aol::ParameterParser pparser( parsername ); 
		
		int numberOfshell = pparser.getInt("numberOfshell");
		
		RType gamma1 = pparser.getDouble("memWeight"); 
		RType gamma2 = pparser.getDouble("bendWeight"); 
		
		string srcDir(pparser.getString("srcDirectory") );
		string destDir( pparser.getString("destDirectory") );
		string filestem( pparser.getString("loadFileStem") );
		string savestem( pparser.getString("saveFileStem") );
		//! load reference mesh
		ostringstream nameOfRefMesh;
		nameOfRefMesh << srcDir << filestem << 1 << ".ply";
		cerr << "Load reference: " << filestem << endl;
		TriMeshType refMesh( nameOfRefMesh.str() );          
		cerr << "No. of Ver " << refMesh.getNumVertices() << endl;
		
		//! init topology saver
		MeshTopologySaver<TriMeshType> Topology( refMesh );
		
		//! load meshes
		aol::VectorContainer< aol::MultiVector<RType>  > geometriesOfInputData( numberOfshell );    
		for( int i = 0; i < numberOfshell; i++ ){
			ostringstream nameOfMesh;
			nameOfMesh << srcDir << filestem << i << ".ply";
			cerr << "Load " << nameOfMesh.str() << endl;
			refMesh.loadFromPLY( nameOfMesh.str() );
			refMesh.toVector( geometriesOfInputData[i] );
		} 
		
		
		//! compute the elastic average S = min ( sum_{i=1} ^ m( s_i - S))
		bool fixBoundary = pparser.checkAndGetBool( "fixBoundary" );
		
		ElasticAverageFunctional< TriMeshType, MembraneDeformation<TriMeshType>, SimpleBendingDeformation<TriMeshType> > elasticEnergy( pparser, Topology,  geometriesOfInputData );  
		ElasticAverageGradient< TriMeshType, MembraneDeformation<TriMeshType>, SimpleBendingDeformation<TriMeshType> > elasticGrad( pparser, Topology,  geometriesOfInputData );      
		ElasticAverageHessian< TriMeshType, MembraneDeformation< TriMeshType>, SimpleBendingDeformation < TriMeshType > > elasticHessian( pparser, Topology,  geometriesOfInputData );
		
		//! fix boundary?
		aol::BitVector BdryMask( Topology.getNumVertices() );
		if( fixBoundary ){
			cerr << "Fix boundary.\n";
			refMesh.fillBoundaryMask( BdryMask );
			elasticGrad.setBoundaryMask( BdryMask );
			elasticHessian.setBoundaryMask( BdryMask );
		}   
		
		//! initialize with linear average
		aol::MultiVector<RType> InitialValue( geometriesOfInputData[0], aol::STRUCT_COPY );
		InitialValue.setZero();
		for( int i = 0; i < numberOfshell; i++ )
			InitialValue.addMultiple( geometriesOfInputData[i], 1. / static_cast<RType>( numberOfshell ) );
		
		aol::MultiVector<RType> ElasticAverage( geometriesOfInputData[0], aol::STRUCT_COPY );
		ElasticAverage = InitialValue;
		//! save init ElasticAverage
		ostringstream nameOfInitAve;
		nameOfInitAve << destDir << "initAve.ply";
		refMesh.fromVector( ElasticAverage );
		refMesh.saveAsPLY(nameOfInitAve.str());
		
		//! Compute the nodal difference of u_i = \tilde {S_i} - S
		for( int i = 0; i < numberOfshell; i++ ){
			geometriesOfInputData[i] -= ElasticAverage;
		}
		
		//! Assemble the Hessian matrix W_S [S], using the assembleHessianDef
		cerr << "Start assembling dot product...\n";
		aol::SparseBlockMatrix< aol:: SparseMatrix<RType> > HessMat ( 3, 3 );
		for ( int i = 0; i < 3; i++ )
			for ( int j = 0; j < 3; j++ )
			{
				HessMat.allocateMatrix ( i, j, refMesh.getNumVertices(), refMesh.getNumVertices() );
			}
			
		MembraneDeformation<TriMeshType> memDeformation( Topology, pparser );
		memDeformation.assembleAddDefHessian( ElasticAverage, ElasticAverage, HessMat, gamma1 );
	
		SimpleBendingDeformation < TriMeshType > bendDeformation( Topology, pparser ); 
		bendDeformation.assembleAddDefHessian( ElasticAverage, ElasticAverage, HessMat, gamma2 );
		
		//! Do PCA using aol::pcaOp ( dotProd, ( u_i ))
		// perform the PCA 
		// caution: covariance operator does not divide by (number of samples - 1) !
		cerr << "Start pca...\n";
		aol::VectorContainer<aol::MultiVector<RType> > modes;
		aol::Vector<RType> variances;
		aol::PCAOp<RType,aol::MultiVector<RType> > pca( HessMat, numberOfshell );
		pca.modesAndVariances( geometriesOfInputData, modes, variances );
		
		//! display variances
		cerr << "Variances = " << variances << endl;
		
		//1 Visualize results
		cerr << "Start storing...\n";
		for( int i = 0; i < numberOfshell; i++ ){      
			for( int j = 0; j < 3; j++ ){
				aol::MultiVector<RType> temp( ElasticAverage );
				temp.addMultiple( modes[i], 1.0 * j );
				refMesh.fromVector( temp );
				ostringstream nameOfModeMesh;
				nameOfModeMesh << destDir << savestem << i << "_" << j << ".ply";
				refMesh.saveAsPLY( nameOfModeMesh.str() );
			}
			
			for( int j = 1; j < 3; j++ ){
				aol::MultiVector<RType> temp( ElasticAverage );
				temp.addMultiple( modes[i], -1.0 * j );
				refMesh.fromVector( temp );
				ostringstream nameOfModeMesh;
				nameOfModeMesh << destDir << savestem << i << "_-" << j << ".ply";
				refMesh.saveAsPLY( nameOfModeMesh.str() );
			}
		}   
	} catch ( aol::Exception &el){
		el.dump();
	}
	
	aol::callSystemPauseIfNecessaryOnPlatform();
	return 0;
} 