//Function checks the neighbours and if it is not in the list adds it

void addNeighbours 
    (
        DynamicList<label>& neigh,
        const label& airCell,
        const Foam::fvMesh& mesh,
        const volScalarField& alpha1,
        const float& alphaLim,
        volScalarField& bubble,
        const int& ii
    )
{
    bubble[airCell] = ii;
    const labelList& neighbours = mesh.cellCells()[airCell]; //create the list of neighbours
    
    forAll(neighbours, i) // check every neighbour of initial cell
    {
            if (alpha1[neighbours[i]] > alphaLim)
                {
                    neigh.append(neighbours[i]);
                    bubble[neighbours[i]] = ii;
                }
    }

    forAll (neigh, i) // check every neighbour of neighbour of initial cell
    {
        const labelList& neighOfNeigh = mesh.cellCells()[neigh[i]];
        
        forAll(neighOfNeigh, j) // check each neighbour
        {
            if (alpha1[neighOfNeigh[j]] > alphaLim && bubble[neighOfNeigh[j]] < ii)
            {
                neigh.append(neighOfNeigh[j]);
                bubble[neighOfNeigh[j]] = ii;
            }
        }
    }

}

void calcPos
    (
        const DynamicList<DynamicList<label>>& Nn,
        const Foam::fvMesh& mesh,
        DynamicList<vector>& partPos
    )
{
    forAll(Nn,i)
    {
        vector pos = Zero;
        int n = 0;
        forAll(Nn[i],j)
        {

            pos = pos + mesh.C()[Nn[i][j]];
            ++n;
        }
        pos = pos / n;
        partPos.append(pos);
    }
}

void calcLimits
    (
        const DynamicList<DynamicList<label>>& Nn,
        const Foam::fvMesh& mesh,
        DynamicList<DynamicList<vector>>& partLimits
    )
{
    forAll(Nn,i)
    {
	DynamicList<vector> tmp;
        vector top = Zero;
	vector bottom = Zero;
	bool s = true;
        forAll(Nn[i],j)
        {
	    if(s == true)
	    {
		top = mesh.C()[Nn[i][j]];
		bottom = mesh.C()[Nn[i][j]];
	    }
	   else
	   {
		if(mesh.C()[Nn[i][j]].component(vector::X) > top.component(vector::X)) top.component(vector::X) = mesh.C()[Nn[i][j]].component(vector::X);
		if(mesh.C()[Nn[i][j]].component(vector::Y) > top.component(vector::Y)) top.component(vector::Y) = mesh.C()[Nn[i][j]].component(vector::Y);
		if(mesh.C()[Nn[i][j]].component(vector::Z) > top.component(vector::Z)) top.component(vector::Z) = mesh.C()[Nn[i][j]].component(vector::Z);

                if(mesh.C()[Nn[i][j]].component(vector::X) < bottom.component(vector::X)) bottom.component(vector::X) = mesh.C()[Nn[i][j]].component(vector::X);
                if(mesh.C()[Nn[i][j]].component(vector::Y) < bottom.component(vector::Y)) bottom.component(vector::Y) = mesh.C()[Nn[i][j]].component(vector::Y);
                if(mesh.C()[Nn[i][j]].component(vector::Z) < bottom.component(vector::Z)) bottom.component(vector::Z) = mesh.C()[Nn[i][j]].component(vector::Z);
	   }
            s = false;
        }
        tmp.append(bottom);
	tmp.append(top);
	partLimits.append(tmp);
    }
}




void calcVolume
    (
        const DynamicList<DynamicList<label>>& Nn,
        const Foam::fvMesh& mesh,
        const volScalarField& alpha1,
        DynamicList<scalar>& partVol
    )
{
    forAll(Nn,i)
    {
        scalar Vol = Zero;
        forAll(Nn[i],j)
        {

            Vol = Vol + (mesh.V()[Nn[i][j]]* (1-alpha1[Nn[i][j]]));
            //Vol = Vol + mesh.V()[Nn[i][j]];
        }
        partVol.append(Vol);
    }
}

void calcDiam
    (
        const DynamicList<scalar>& partVol,
        DynamicList<scalar>& partDiam
    )
{
    forAll(partVol,i)
    {
        
        scalar Diam = std::pow(std::abs(1.90985935*partVol[i]),0.333333);// cbrt();
        partDiam.append(Diam);
    }
}

void calcVel
    (
        const DynamicList<DynamicList<label>>& Nn,
        const Foam::fvMesh& mesh,
        const volVectorField& U,
        DynamicList<vector>& partVel
    )
{
    forAll(Nn,i)
    {
        vector Vel = Zero;
        int n = 0;
        forAll(Nn[i],j)
        {

            Vel = Vel + U[Nn[i][j]];
            ++n;
        }
        Vel = (Vel / n);
        partVel.append(Vel); //<-just Vel
    }
}

void printVector
    (
        const vector& vec,
        std::ofstream& file,
        const string& varName
    )
{
    file << varName << '\t' << "(" << vec.component(vector::X) << " " << vec.component(vector::Y) << " " << vec.component(vector::Z) << ")" << ";" << std::endl << std::endl;
}

void printVector
    (
        const vector& vec,
        std::ofstream& file
    )
{
    file  << "(" << vec.component(vector::X) << " " << vec.component(vector::Y) << " " << vec.component(vector::Z) << ")";
}

void printFakeVector
    (
        const vector& vec,
        const scalar& Xpos,
        std::ofstream& file
    )
{
    file  << "(" << vec.component(vector::X)  << " " << vec.component(vector::Y) << " " << Xpos << ")";
}


void printVectorList
    (
        const DynamicList<vector>& vecList,
        std::ofstream& file,
        const string& varName,
        const DynamicList<label> activeList
    )
{
    forAll(activeList,i)
    {
        file << varName << i << '\t';
        printVector(vecList[activeList[i]], file);
        file << ";" <<std::endl;
    }
    file << std::endl;
}

void printFakePosVectorList
    (
        const DynamicList<vector>& vecList,
        std::ofstream& file,
        const string& varName,
        const DynamicList<label> activeList
    )
{
    forAll(activeList,i)
    {
        file << varName << i << '\t';
        printFakeVector(vecList[activeList[i]], 0.00025 ,file);
        file << ";" <<std::endl;
    }
    file << std::endl;
}

void printScalarList
    (
        const DynamicList<scalar>& scalarList,
        std::ofstream& file,
        const string& varName,
        const DynamicList<label> activeList
    )
{
    forAll(activeList,i)
    {
        file << varName << i << '\t' << scalarList[activeList[i]]<< ";" << std::endl;
    }
    file << std::endl;
}

DynamicList<label> activeParticles
(
    const DynamicList<vector>& partPos
)
{
    DynamicList<label> list;
    //filter all particles after filter area
    forAll(partPos, i)
    {
        list.append(i);
    }
    return list;
}

