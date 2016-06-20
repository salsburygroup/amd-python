from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship
from .Database import Base

class Docking(Base):
    __tablename__ = 'docking'
    id = Column(Integer, primary_key=True)
    ligand_name = Column(String(250), nullable=False)
    receptor_name = Column(String(250), nullable=False)

    def __init__(self, ligand_name, receptor_name):
        self.ligand_name = ligand_name
        self.receptor_name = receptor_name

    def __repr__(self):
        return 'ligand: {0} receptor: {1}'.format(self.ligand_name, self.receptor_name)


class Energy(Base):
    __tablename__ = 'energy'
    id = Column(Integer, primary_key=True)
    rank = Column(Integer, nullable=False)
    dG = Column(Float, nullable=False)
    docking_id = Column(Integer, ForeignKey('docking.id'))
    docking = relationship(Docking)
