/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *  Side Effects Software Inc
 *  123 Front Street West
 *  Toronto, Ontario
 *  Canada   M5V 3E7
 *  416-504-9876
 *
 * NAME:    crowd_relationships.h
 */

#ifndef __crowd_relationships_h__
#define __crowd_relationships_h__

#include <math.h>

#define CROWD_REL_CONSTRAIN_T 1
#define CROWD_REL_CONSTRAIN_R 2
#define CROWD_REL_CONSTRAIN_ALL (CROWD_REL_CONSTRAIN_T | CROWD_REL_CONSTRAIN_R)

struct AgentXform
{
    vector myT;
    vector4 myR;
    vector myS;

    void set_identity()
    {
        myT = {0, 0, 0};
        myR = {0, 0, 0, 1};
        myS = {1, 1, 1};
    }

    void set_xform(const vector t; const vector4 r; const vector s)
    {
        myT = t;
        myR = r;
        myS = s;
    }

    void set_xform(const vector t; const vector4 r)
    {
        this->set_xform(t, r, {1, 1, 1});
    }

    void set_from_matrix(const matrix xform)
    {
        vector r;
        cracktransform(XFORM_SRT, XFORM_XYZ, {0, 0, 0}, xform, myT, r, myS);
        myR = eulertoquaternion(radians(r), XFORM_XYZ);
    }

    void get_xform(export vector t; export vector4 r; export vector s)
    {
        t = myT;
        r = myR;
        s = myS;
    }

    void mult(const AgentXform xform)
    {
        myT *= xform.myS;
        myT = qrotate(xform.myR, myT) + xform.myT;
        myR = qmultiply(xform.myR, myR);
    }

    AgentXform compute_difference(const AgentXform xform)
    {
        vector4 delta_r = qmultiply(qinvert(myR), xform.myR);
        vector delta_t = xform.myT - qrotate(delta_r, myT);

        AgentXform delta;
        delta->set_xform(delta_t, delta_r);
        return delta;
    }

    void print()
    {
        printf("t: %g, r: %g, s: %g\n", myT, myR, myS);
    }
};

struct AgentNode
{
    void init(const int ptnum; const int parent_idx)
    {
        myPtNum = ptnum;
        myPrimNum = pointprims(0, ptnum)[0];
        myParentIdx = parent_idx;

        AgentXform xform;
        xform->set_from_matrix(
            primintrinsic(0, "packedfulltransform", myPrimNum));
        xform->get_xform(myT, myR, myS);
    }

    AgentXform get_agent_xform()
    {
        AgentXform xform;
        xform->set_xform(myT, myR, myS);
        return xform;
    }

    AgentXform
    build_constraint_xform(const AgentNode node; const string joint_name;
                           const int constraint_type)
    {
        int joint_idx = agentrigfind(0, node.myPrimNum, joint_name);

        AgentXform xform;
        if (joint_idx < 0)
            xform->set_identity();
        else
        {
            xform->set_from_matrix(
                agentworldtransform(0, node.myPrimNum, joint_idx));

            if (!(constraint_type & CROWD_REL_CONSTRAIN_T))
                xform.myT = {0, 0, 0};

            if (!(constraint_type & CROWD_REL_CONSTRAIN_R))
                xform.myR = {0, 0, 0, 1};
        }

        xform->mult(node->get_agent_xform());

        return xform;
    }

    void transform(const AgentNode parent)
    {
        int constraint_type = point(0, "agentrel_type", this.myPtNum);
        if (!constraint_type)
            constraint_type = CROWD_REL_CONSTRAIN_ALL;

        string parent_joint_name = point(0, "agentrel_parentjoint", myPtNum);
        AgentXform parent_xform = this->build_constraint_xform(
            parent, parent_joint_name, constraint_type);

        string child_joint_name = point(0, "agentrel_childjoint", myPtNum);
        AgentXform child_xform = this->build_constraint_xform(
            this, child_joint_name, constraint_type);

        // Determine the change needed to align the child with the parent.
        AgentXform align_with_parent =
            child_xform->compute_difference(parent_xform);

        // Adjust the agent's overall transform.
        AgentXform agent_xform = this->get_agent_xform();
        agent_xform->mult(align_with_parent);
        agent_xform->get_xform(myT, myR, myS);

        // Apply the extra offset.
        {
            vector offset_t = point(0, "agentrel_P", this.myPtNum);

            int success;
            vector4 offset_r =
                pointattrib(0, "agentrel_orient", this.myPtNum, success);
            if (!success)
                offset_r = {0, 0, 0, 1};

            myT += qrotate(myR, offset_t);
            myR = qmultiply(myR, offset_r);
        }
    }

    void save()
    {
        setpointattrib(0, "P", myPtNum, myT);
        addpointattrib(0, "orient", {0, 0, 0, 1}, "quaternion");
        setpointattrib(0, "orient", myPtNum, myR);
        setpointattrib(0, "pscale", myPtNum, max(myS));

        // If pointinstancetransform is enabled (which it is by default for
        // agents), clear the primitive transform since we're describing it
        // with the standard point attributes. Otherwise, we need to update the
        // transform intrinsic (this lets us also attach normal packed
        // primitives, which might be useful).
        matrix3 prim_xform = ident();
        if (!primintrinsic(0, "pointinstancetransform", myPrimNum))
        {
            scale(prim_xform, myS);
            prim_xform *= qconvert(myR);
        }

        setprimintrinsic(0, "transform", myPrimNum, prim_xform);
    }

    int id()
    {
        int success;
        int agent_id = pointattrib(0, "id", myPtNum, success);
        return success ? agent_id : myPtNum;
    }

    int myPtNum;
    int myPrimNum;
    int myParentIdx; // Index in the hierarchy array.

    // Note: for now, we can't store an AgentXform struct here due to a crash
    // in VEX.
    vector myT;
    vector4 myR;
    vector myS;
};

/// Update the transforms of this agent's children.
void
crowdrelationship_updatetransforms(const int root_ptnum)
{
    AgentNode root_agent;
    root_agent->init(root_ptnum, -1);

    AgentNode hierarchy[];
    push(hierarchy, root_agent);

    for (int i = 0; i < len(hierarchy); ++i)
    {
        AgentNode agent = hierarchy[i];

        // Apply the parent's transform.
        if (agent.myParentIdx >= 0)
        {
            AgentNode parent_agent = hierarchy[agent.myParentIdx];
            agent->transform(parent_agent);
            hierarchy[i] = agent;
        }

        // Add children.
        int agent_id = agent->id();
        if (agent_id < 0)
            continue;

        int n = findattribvalcount(0, "point", "agentrel_parentid", agent_id);
        for (int j = 0; j < n; ++j)
        {
            int child_ptnum =
                findattribval(0, "point", "agentrel_parentid", agent_id, j);

            AgentNode child_xform;
            child_xform->init(child_ptnum, i);
            push(hierarchy, child_xform);
        }
    }

    // Write out the updated transforms.
    foreach (AgentNode agent; hierarchy)
    {
        agent->save();
    }
}

#endif
