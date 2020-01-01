/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::incompressible::LESModels::SchumannGrotzbach_TimeAvgFvPatchField

Description:
    Wall boundary condition for total shear stress at a solid boundary. Stresses
    are computed using the Schumann/Grotzbach formulation. In any inhomogeneous
    surface, average velocity is obtained though a time average instead of a
    horizontal-average or local mode.

SourceFiles
    SchumannGrotzbach_TimeAvgFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef SchumannGrotzbach_TimeAvgFvPatchField_H
#define SchumannGrotzbach_TimeAvgFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class SchumannGrotzbach_TimeAvgFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class SchumannGrotzbach_TimeAvgFvPatchField
:
    public fixedValueFvPatchSymmTensorField
{

private:

    // Private data

        //- Von Karman constant
        scalar kappa_;

        //- Roughness height
        scalarField z0_;

        //- Constant for unstable Monin-Obuhkov scaling law
        scalar betaM_;

        // -Constant for stable Monin-Obuhkov scaling law
        scalar gammaM_;

        // -Type of application of Monin-Obuhkov scaling, local or nonlocal
        word averageType_;

	//- Surface Velocity Field // Add by Y.Han================================
        vectorField  U1_Paral_local_timeAve_;

	//
	bool restart_;

	//- Running time-average range

	scalar TimeRange_;

	//- Current starting index in sampleTimes
        label startSampleTime_;
        scalar StartTime_;

    // Private member functions

        // -Function to transform vectors from Cartesian into terrain-local
        //  coordinates
        vector transformVectorCartToLocal
        (
             vector v,
             vector xP,
             vector yP,
             vector zP
        );

        // -Function to transform symmetric tensors from terrain-local to
        //  Cartesian coordinates
        symmTensor transformSymmTensorLocalToCart
        (
             symmTensor SP,
             vector xP,
             vector yP,
             vector zP
        );

        // -Function to calculate friction velocity using Monin-Obuhkov
        //  similarity laws
        void uStarEvaluate
        (
             scalar &uStar,
             scalar& L,
             scalar& phiM,
             scalar U,
             scalar z1,
             scalar z0,
             scalar kappa,
             scalar gammaM,
             scalar betaM,
             scalar g,
             scalar TRef,
             scalar qw,
             scalar eps,
             scalar tol,
             label iterMax
        );


public:

    //- Runtime type information
    TypeName("SchumannGrotzbach_TimeAvg");


    // Constructors

        //- Construct from patch and internal field
        SchumannGrotzbach_TimeAvgFvPatchField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SchumannGrotzbach_TimeAvgFvPatchField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  stressABLWallFunctionFvPatchScalarField
        //  onto a new patch
        SchumannGrotzbach_TimeAvgFvPatchField
        (
            const SchumannGrotzbach_TimeAvgFvPatchField&,
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        SchumannGrotzbach_TimeAvgFvPatchField
        (
            const SchumannGrotzbach_TimeAvgFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchSymmTensorField> clone() const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new SchumannGrotzbach_TimeAvgFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SchumannGrotzbach_TimeAvgFvPatchField
        (
            const SchumannGrotzbach_TimeAvgFvPatchField&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchSymmTensorField> clone
        (
            const DimensionedField<symmTensor, volMesh>& iF
        ) const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new SchumannGrotzbach_TimeAvgFvPatchField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

    	// Mapping functions

	//- Reverse map the given fvPatchField onto this fvPatchField
	virtual void rmap
	(
	    const fvPatchSymmTensorField&,
	    const labelList&
	);

        // Access functions
            //- Return z_0
            scalarField& z0()
            {
                return z0_;
            }

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
