IOdictionary myDict
(

    IOobject
    (
        "constant/myDictFile",
        "",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )

);

