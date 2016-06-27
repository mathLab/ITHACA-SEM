<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fourier Single Mode Adjoint Basis, P=7, Baseflow file</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>SM_Adj_base_file.xml</parameters>
    <files>
        <file description="Session File">SM_Adj_base_file.xml</file>
        <file description="Base Flow File">SM_Adj_base_file.bse</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.13451e-12</value>
            <value variable="v" tolerance="1e-12">1.34118e-13</value>
            <value variable="w" tolerance="1e-12">6.3238e-13</value>
            <value variable="p" tolerance="1e-12">4.76387e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.85006e-12</value>
            <value variable="v" tolerance="1e-12">8.81739e-13</value>
            <value variable="w" tolerance="1e-12">1.78049e-12</value>
            <value variable="p" tolerance="1e-12">1.43926e-10</value>
        </metric>
    </metrics>
</test>
