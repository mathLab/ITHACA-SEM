<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=4, periodic BCs</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_m4_per_hdf5.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Tet_channel_m4_per_hdf5.xml</file>
        <file description="Session File">Tet_channel_m4_per_hdf5.nekg</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">8.59656e-14</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-8">2.50344e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">7.42739e-14</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-8">1.50102e-13</value>
        </metric>
    </metrics>
</test>
