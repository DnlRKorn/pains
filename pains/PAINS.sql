CREATE TABLE pains_compounds(
CID integer PRIMARY KEY,
SMILES varchar(255),
PAINS_A varchar(31),
PAINS_B varchar(31),
PAINS_C varchar(31),
AllAssays_Active integer,
AllAssays_Total integer,
LuciferaseAssays_Active integer,
LuciferaseAssays_Total integer,
BetaLactamaseAssays_Active integer,
BetaLactamaseAssays_Total integer,
FluorescenceAssays_Active integer,
FluorescenceAssays_Total integer
);

