//
// Created by Rina Karnauch & Ofek Kaveh on 14/04/2021.
//


/**
 * main for the struct alignment
 * @param argc amount of args to the program
 * @param argv array of arguments as chars
 * @return none
 */
#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <iostream>
#include "Triangle.h"

using namespace std;


vector<Triangle> makeCombinations(const Molecule<Atom> &molecule)
/**
 * method to get all triangle combinations out of atom
 * @param molecule: current molecule to turn into triangles
 * @return vector of all triangles of conseccutive CA triangles
 */
{
	vector<Triangle> triangles;
	for (int i = 0; i < molecule.size() - 2; i++)
	{

		Atom x = molecule[i];
		Atom y = molecule[i + 1];
		Atom z = molecule[i + 2];
		Triangle t = Triangle(x, y, z);
		triangles.push_back(t);
	}
	return triangles;
}

RigidTrans3 getTransformation(const Triangle& model, const Triangle& target)
/**
 * method to get transformation between two triangles
 * @param model model triangle to put onto
 * @param target target triangle to put model ontop of it
 * @return a rigidtrans of model->target
 */
{
	return target | model;
}

vector<RigidTrans3>
getAllTransformations(vector<Triangle> mulModel,
					  vector<Triangle> mulTarget)
/**
 * method to get all possible transformations out of all triangles
 * @param mulModel vector of all model triangles
 * @param mulTarget vector of all target triangles
 * @return vector of all transformation
 */
{
	vector<RigidTrans3> transformations;
	for (int i = 0; i < mulModel.size(); i++)
	{
		for (int j = 0; j < mulTarget.size(); j++)
		{
			transformations.push_back(getTransformation(mulModel[i], mulTarget[j]));
		}
	}
	return transformations;

}



int main(int argc, char *argv[])
/**
 * main to apply transformations and find best
 * @return 0 for success 1 for failure
 */
{
    // measure the run time
//    auto start = std::chrono::system_clock::now();

    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " dist_threshold pdb1 pdb2"
                  << std::endl;
        exit(1);
    }

    // ********************Parameters********************
    float m_fDistThr = atof(argv[1]); // distance threshold on atoms in correspondence

//    std::cout << "Distance threshold: " << m_fDistThr << std::endl;

    // read the two files into Molecule
    Molecule<Atom> molModel, molTarget, molModelFull;

    std::ifstream fileModel(argv[3]);
    std::ifstream fileTarget(argv[2]);
    std::ifstream savedModel(argv[3]);
    std::ifstream tempRead(argv[3]);

    if (!fileModel)
    {
        std::cout << "File " << argv[3] << "does not exist." << std::endl;
        return 0;
    }
    if (!fileTarget)
    {
        std::cout << "File " << argv[2] << "does not exist." << std::endl;
        return 0;
    }

    std::string lineRead;
    std::getline(tempRead, lineRead);
    savedModel.clear();
    molModelFull.readPDBfile(savedModel, PDB::AllSelector());

	// if first line is RNA
    if (lineRead.find("RNA") != std::string::npos)
    {
        molModel.readPDBfile(fileModel, PDB::PSelector());
        molTarget.readPDBfile(fileTarget, PDB::PSelector());
    }
    else
    {
    	// otherwise we read it as a protein
        molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
        molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
    }


    // calculate center of mass
    Vector3 vectModelMass(0, 0, 0);
    for (unsigned int i = 0; i < molModel.size(); i++)
    {
        vectModelMass += molModel[i].position();
    }
    vectModelMass /= molModel.size();

    Vector3 vectTargetMass(0, 0, 0);
    for (unsigned int i = 0; i < molTarget.size(); i++)
    {
        vectTargetMass += molTarget[i].position();
    }
    vectTargetMass /= molTarget.size();

    // transform the molecules to the center of the coordinate system
    molModel += (-vectModelMass);
    molTarget += (-vectTargetMass);

	// next we insert the target molecule into hash
	// this will help us to find atoms that are close faster
	GeomHash<Vector3, int> gHash(3, m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
	for (unsigned int i = 0; i < molTarget.size(); i++)
	{
		gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
	}

    // now we try rotations and choose the best alignment from random rotations
    unsigned int iMaxSize = 0;
    RigidTrans3 rtransBest;
    float RMSD;

    vector<Triangle> modelTriangles = makeCombinations(molModel);
    vector<Triangle> targetTriangles = makeCombinations(molTarget);

    vector<RigidTrans3> transformations = getAllTransformations(modelTriangles, targetTriangles);

    for (int k = 0; k < transformations.size(); k++)
    {
        Match match;
        // apply rotation on each atom in the model molecule and
        // add the pairs of atoms (one from target and one from model)
        // that are close enough to the match list
        for (unsigned int i = 0; i < molModel.size(); i++)
        {
            Vector3 mol_atom = transformations[k] * molModel[i].position(); // rotate

            // find close target molecule atoms using the hash
            HashResult<int> result;
            gHash.query(mol_atom, m_fDistThr, result); // key is mol atom coordinate

            // check if the atoms in the result are inside the distance threshold
            // the hash is a cube shape, there can be atoms further that the threshold
            for (auto x = result.begin(); x != result.end(); x++)
            {
                float dist = mol_atom.dist(molTarget[*x].position());
                if (dist <= m_fDistThr)
                {
                    float score = (1 / (1 + dist));
                    match.add(*x, i, score, score);
                }
            }
            result.clear();
        }

        match.calculateBestFit(molTarget, molModel);

        if (iMaxSize < match.size())
        {
            iMaxSize = match.size();
            rtransBest = match.rigidTrans();
            RMSD = match.rmsd();
        }
    }

    std::cout << iMaxSize << " " << RMSD << " " << rtransBest << " " << std::endl;


	molModelFull *=rtransBest;

	std::ofstream fileStream("transformed.pdb", std::ofstream::out);
	fileStream << molModelFull;
}



