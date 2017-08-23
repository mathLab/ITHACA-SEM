<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, axi-symmetric nozzle without swirl</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Nozzle_AxiSym_NoSwirl.xml</parameters>
    <files>
        <file description="Session File">Nozzle_AxiSym_NoSwirl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.16488</value>
            <value variable="rhou" tolerance="1e-12">0.593021</value>
            <value variable="rhov" tolerance="1e-12">43.6473</value>
            <value variable="E" tolerance="1e-12">643932</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.25895</value>
            <value variable="rhou" tolerance="1e-12">0.732933</value>
            <value variable="rhov" tolerance="1e-12">44.9363</value>
            <value variable="E" tolerance="1e-12">260029</value>
        </metric>
    </metrics>
</test>


