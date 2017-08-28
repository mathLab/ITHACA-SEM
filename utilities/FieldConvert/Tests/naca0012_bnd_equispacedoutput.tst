<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Equispaced output of a 2DH1D field </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m equispacedoutput naca0012_bnd.xml naca0012_b0.fld output.vtu</parameters>
    <files>
        <file description="Session File">naca0012_bnd.xml</file>
        <file description="Session File">naca0012_b0.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.479474</value>
            <value variable="y" tolerance="1e-6">0.0410293</value>
            <value variable="z" tolerance="1e-6">0.586302</value>
            <value variable="u" tolerance="1e-12">3.92345e-13</value>
            <value variable="v" tolerance="1e-12">9.11782e-13</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0.330907</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">0.0600037</value>
            <value variable="z" tolerance="1e-6">1</value>
            <value variable="u" tolerance="1e-12">1.22759e-12</value>
            <value variable="v" tolerance="1e-12">7.24523e-13</value>
            <value variable="w" tolerance="1e-12">1.06335e-24</value>
            <value variable="p" tolerance="1e-12">0.685352</value>
        </metric>
    </metrics>
</test>

