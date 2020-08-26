<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Structured grid generation from NekMesh wrapper</description>
    <executable python="true">LoadCAD.py</executable>
    <parameters>3d_sphere.stp output.xml</parameters>
    <files>
        <file description="STEP file input">../../../../../utilities/NekMesh/Tests/MeshGen/STEP/3d_sphere.stp</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
