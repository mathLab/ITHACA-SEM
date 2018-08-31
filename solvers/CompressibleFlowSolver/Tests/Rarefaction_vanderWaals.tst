<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Rarefaction wave, van der Waals equation of state</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Rarefaction_vanderWaals.xml</parameters>
    <files>
        <file description="Session File">Rarefaction_vanderWaals.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">189.173</value>
            <value variable="rhou" tolerance="1e-12">12445.1</value>
            <value variable="rhov" tolerance="1e-12">36.3746</value>
            <value variable="E" tolerance="1e-12">2.06388e+08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">202.888</value>
            <value variable="rhou" tolerance="1e-12">13348</value>
            <value variable="rhov" tolerance="1e-12">1460.79</value>
            <value variable="E" tolerance="1e-12">2.21351e+08</value>
        </metric>
    </metrics>
</test>

