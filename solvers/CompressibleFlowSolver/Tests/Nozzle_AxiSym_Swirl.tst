<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, axi-symmetric nozzle with swirl</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Nozzle_AxiSym_Swirl.xml</parameters>
    <files>
        <file description="Session File">Nozzle_AxiSym_Swirl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.16478</value>
            <value variable="rhou" tolerance="1e-12">0.591809</value>
            <value variable="rhov" tolerance="1e-12">43.6375</value>
            <value variable="rhow" tolerance="1e-12">4.69675</value>
            <value variable="E" tolerance="1e-12">643926</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.25873</value>
            <value variable="rhou" tolerance="1e-12">0.732933</value>
            <value variable="rhov" tolerance="1e-12">44.9363</value>
            <value variable="rhow" tolerance="1e-12">14.2803</value>
            <value variable="E" tolerance="1e-12">260113</value>
        </metric>
    </metrics>
</test>


