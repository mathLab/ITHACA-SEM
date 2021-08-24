<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Flow forced by a derivative discontinuity along a line with coupled solver, P=5</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>2DFlow_lineforcing_bcfromfile.xml</parameters>
    <files>
        <file description="Session File">2DFlow_lineforcing_bcfromfile.xml</file>
        <file description="Session File">2DFlow_lineforcing_bcfromfile_u_1.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.612387</value>
            <value variable="v" tolerance="1e-6">2.96601e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">1.5</value>
            <value variable="v" tolerance="1e-6">1.25209e-13</value>
        </metric>
    </metrics>
</test>


