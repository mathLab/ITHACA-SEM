<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, LE simulation, BCs from file, FRHU, SEM, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>RAE5240_BSF_LE_bcsfromfile_par.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">RAE5240_BSF_LE_bcsfromfile_par.xml</file>
        <file description="Session File">RAE5240_BSF_LE_bcsfromfile_par.rst</file>
        <file description="Session File">RAE5240_BSF_LE_bcsfromfile_par_bc1.bc</file>
        <file description="Session File">RAE5240_BSF_LE_bcsfromfile_par_bc2.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.5526</value>
            <value variable="rhou" tolerance="1e-12">348.093</value>
            <value variable="rhov" tolerance="1e-8">21.5936</value>
            <value variable="E" tolerance="1e-12">358698</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.52341</value>
            <value variable="rhou" tolerance="1e-12">298.756</value>
            <value variable="rhov" tolerance="1e-8">164.223</value>
            <value variable="E" tolerance="1e-12">349234</value>
        </metric>
    </metrics>
</test>


