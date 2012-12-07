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
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.58914</value>
            <value variable="v" tolerance="1e-12">0.00237598</value>
            <value variable="p" tolerance="1e-12">0.00236396</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1</value>
            <value variable="v" tolerance="1e-12">0.00200964</value>
            <value variable="p" tolerance="1e-12">0.00177702</value>
        </metric>
    </metrics>
</test>


