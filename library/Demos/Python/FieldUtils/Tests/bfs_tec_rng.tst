<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D tecplot output with a range restriction</description>
    <executable python="true">bfs_tec_rng.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1.35354</value>
            <value variable="y" tolerance="1e-6">1.09498</value>
            <value variable="u" tolerance="1e-6">1.15931</value>
            <value variable="v" tolerance="1e-6">0.0195245</value>
            <value variable="p" tolerance="1e-6">0.119152</value>
        </metric>
    </metrics>
</test>

