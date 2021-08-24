<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Adjoint stability (Arpack): Channel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanStability_adj_Ar.xml</parameters>
    <files>
        <file description="Session File">ChanStability_adj_Ar.xml</file>
	<file description="Session File">ChanStability_adj_Ar.bse</file>
	<file description="Session File">ChanStability_adj_Ar.rst</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.00031,0.0349782</value>
            <value tolerance="0.001">1.00031,-0.0349782</value>
        </metric>
    </metrics>
</test>


