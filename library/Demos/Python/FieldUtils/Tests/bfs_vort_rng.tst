<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D vorticity output </description>
    <executable python="true">bfs_vort_rng.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-3">1.35354</value>
            <value variable="y" tolerance="1e-3">1.09498</value>
            <value variable="u" tolerance="1e-3">1.15931</value>
            <value variable="v" tolerance="1e-3">0.019524</value>
            <value variable="p" tolerance="1e-3">0.119152</value>
            <value variable="W_z" tolerance="1e-3">3.47993</value>
        </metric>
    </metrics>
</test>
